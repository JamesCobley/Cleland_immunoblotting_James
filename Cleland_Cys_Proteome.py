import pandas as pd

# Load the Excel file
file_path = '/content/Supplementary data file 5.xlsx'  # Replace with the correct path
xls = pd.ExcelFile(file_path)

# List to store measurable protein data
measurable_proteins = []

# Loop through each sheet in the Excel file
for sheet in xls.sheet_names:
    df = pd.read_excel(xls, sheet)
    
    # Ignore columns A-C which correspond to 0 cysteine residues
    # We will loop over the cysteine columns starting from where cysteine counts begin
    for col in df.columns:
        if 'Cysteine_Residue_Count' in col and '0_' not in col:  # Skip if cysteine count is 0
            cysteine_count_col = col
            cysteine_mass_col = col.replace('Residue_Count', 'Molecular_Mass_kDa')
            cysteine_id_col = col.replace('Residue_Count', 'IDs')  # Assuming UniProt ID column
            
            if cysteine_mass_col in df.columns and cysteine_id_col in df.columns:
                # Iterate through each protein entry
                for i, row in df.iterrows():
                    cysteine_count = row[cysteine_count_col]
                    molecular_mass = row[cysteine_mass_col] / 1000  # Convert to kDa from Daltons
                    uniprot_id = row[cysteine_id_col]  # UniProt ID
                    
                    # Calculate the oxidized mass
                    oxidized_mass = molecular_mass + (cysteine_count * 5)
                    
                    # Check if protein is measurable
                    if molecular_mass <= 150 and oxidized_mass <= 200:
                        # Add to measurable proteins list
                        measurable_proteins.append({
                            'UniProt_ID': uniprot_id,
                            'Cysteine_Residue_Count': cysteine_count,
                            '100%-Reduced_Molecular_Mass': molecular_mass,
                            '100%-Oxidised_Molecular_Mass': oxidized_mass
                        })

# Create a DataFrame for measurable proteins
measurable_df = pd.DataFrame(measurable_proteins)

# Write to an Excel file
output_file = '/content/measurable_proteins.xlsx'  # Output file path
measurable_df.to_excel(output_file, index=False)

# Print the number of measurable proteins
print(f"Number of measurable proteins: {len(measurable_proteins)}")
print(f"Excel file with measurable proteins saved to: {output_file}")
