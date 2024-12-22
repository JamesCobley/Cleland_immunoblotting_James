import numpy as np
import matplotlib.pyplot as plt

# Define k-space values and their contributions
k_values = np.array([0, 1, 2, 3])
contributions = np.array([0.12, 0.38, 0.38, 0.12])  # Equal contribution for simplicity

# Generate x values for the wave
x = np.linspace(-1, 4, 500)

# Create the convoluted wave as a sum of all contributions
convoluted_wave = sum(contribution * np.exp(-((x - k) ** 2) / (2 * 0.1 ** 2)) for k, contribution in zip(k_values, contributions))

# Create the bar chart data (deconvoluted components)
deconvoluted_contributions = contributions * 100  # Convert to percentages

# Plot the convoluted wave and deconvoluted bar chart
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Plot 1: Bar chart showing oxidized and reduced percentages
axes[0].bar(["Reduced", "Oxidized"], [50, 50], color=["grey", "black"], alpha=0.7)
axes[0].set_title("Oxidized vs Reduced", fontsize=14)
axes[0].set_ylabel("Percentage (%)", fontsize=12)
axes[0].grid(True, linestyle="--", alpha=0.5)

# Plot 2: Deconvoluted data as a single wave
axes[1].plot(x, convoluted_wave, label="Convoluted 50%-Oxidized Wave", color="black", linewidth=2)
axes[1].set_title("Convoluted Data (50%-Oxidized)", fontsize=14)
axes[1].set_xlabel("k-Space", fontsize=12)
axes[1].set_ylabel("Amplitude", fontsize=12)
axes[1].grid(True, linestyle="--", alpha=0.5)
axes[1].legend(fontsize=10)

# Plot 3: Deconvoluted data as a bar chart
axes[2].bar(k_values, deconvoluted_contributions, color=["black", "black", "black", "black"], alpha=0.7)
axes[2].set_title("Deconvoluted k-Space Contributions", fontsize=14)
axes[2].set_xlabel("k-Space", fontsize=12)
axes[2].set_ylabel("Contribution (%)", fontsize=12)
axes[2].set_xticks(k_values)
axes[2].grid(True, linestyle="--", alpha=0.5)

# Adjust layout and save the figure
plt.tight_layout()
file_path = "/content/convoluted_deconvoluted_oxidized_reduced_final.png"
plt.savefig(file_path, dpi=300)
plt.show()
