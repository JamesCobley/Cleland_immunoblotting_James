import pandas as pd
import matplotlib.pyplot as plt

# Load the original file (Supplementary data file 5) and filtered file (Supplementary data file 6)
original_file = '/content/Supplementary data file 5.xlsx'
filtered_file = '/content/Supplementary data file 6.xlsx'

# Read both files
original_xls = pd.ExcelFile(original_file)
filtered_df = pd.read_excel(filtered_file)

# Dictionary to store counts for each cysteine class
original_counts = {}
filtered_counts = {}

# Loop through each sheet in the original file to count total proteins for each cysteine class
for sheet in original_xls.sheet_names:
    df = pd.read_excel(original_xls, sheet)
    
    # Iterate through the columns to count proteins in each cysteine residue class
    for col in df.columns:
        if 'Cysteine_Residue_Count' in col and '0_' not in col:  # Skip 0 cysteine class
            cysteine_count = col.split('_')[0]  # Extract the cysteine count from column name
            if cysteine_count not in original_counts:
                original_counts[cysteine_count] = 0
            original_counts[cysteine_count] += df[col].notna().sum()  # Count non-null proteins

# Count the proteins in the filtered (measurable) file based on cysteine residue counts
for index, row in filtered_df.iterrows():
    cysteine_count = str(int(row['Cysteine_Residue_Count']))  # Convert to string to match key in original_counts
    if cysteine_count not in filtered_counts:
        filtered_counts[cysteine_count] = 0
    filtered_counts[cysteine_count] += 1

# Calculate the percentage of measurable proteins for each cysteine class
percentages = {}
for cysteine_count in original_counts:
    if cysteine_count in filtered_counts:
        percentages[cysteine_count] = (filtered_counts[cysteine_count] / original_counts[cysteine_count]) * 100
    else:
        percentages[cysteine_count] = 0  # If no measurable proteins in that class

# Prepare data for plotting, grouping cysteine counts >= 50 into "50+"
cysteine_classes = []
percentage_values = []

# Sum all proteins with 50 or more cysteines into the "50+" category
above_50_count = 0
above_50_percentage = 0

for key in sorted([int(k) for k in percentages.keys()]):  # Sort cysteine classes numerically
    if key < 50:
        cysteine_classes.append(str(key))  # Convert all to string to keep types consistent
        percentage_values.append(percentages[str(key)])
    else:
        above_50_count += 1  # Increment for all cysteine classes >= 50
        above_50_percentage += percentages[str(key)]  # Sum the percentages

# Add the "50+" category
if above_50_count > 0:
    cysteine_classes.append("50+")
    percentage_values.append(above_50_percentage / above_50_count)  # Average percentage for the "50+" group

# Plotting the improved bar plot
plt.figure(figsize=(12, 6))

# Create the bar plot
bars = plt.bar(cysteine_classes, percentage_values, color='skyblue')

# Add labels and title with increased font size
plt.xlabel('Cysteine Residue Count', fontsize=14)
plt.ylabel('Percentage of Proteins Amenable', fontsize=14)
plt.title('Percentage of Proteins Amenable by Cysteine Residue Count', fontsize=16)

# Rotate x-axis labels for better readability
plt.xticks(cysteine_classes, rotation=45, fontsize=10)

# Set limits for y-axis (0 to 100%)
plt.ylim(0, 100)

# Add grid for y-axis for better readability
plt.grid(axis='y')

# Add value labels on top of each bar for better clarity
for bar in bars:
    yval = bar.get_height()
    if yval > 0:  # Only label bars with a non-zero height
        plt.text(bar.get_x() + bar.get_width()/2, yval + 1, f'{yval:.1f}%', ha='center', va='bottom', fontsize=10)

# Save the plot as a PNG file at 300 DPI
output_image_path = '/content/percentage_proteins_amenable.png'
plt.tight_layout()
plt.savefig(output_image_path, dpi=300)

# Show the plot
plt.show()

print(f"Plot saved as: {output_image_path}")
