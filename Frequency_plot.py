import numpy as np
import matplotlib.pyplot as plt

# Define binary strings for R = 3 series
binary_strings = ["000", "100", "010", "001", "110", "101", "011", "111"]

# Define x positions for cysteine residues
positions = np.linspace(1, 3, 500)  # Continuous range for smooth waves

# Create a figure for subplots
fig, axes = plt.subplots(2, 4, figsize=(16, 8))
axes = axes.flatten()

# Loop through each binary string and plot its wave
for i, binary in enumerate(binary_strings):
    # Convert binary string to numeric values
    oxidation_states = np.array([int(bit) for bit in binary])
    
    # Create the wave function
    wave = np.zeros_like(positions)
    for j, state in enumerate(oxidation_states):
        if state == 1:
            wave += state * np.exp(-((positions - (j + 1)) ** 2) / 0.1)  # Gaussian peaks
    
    # Plot the wave
    ax = axes[i]
    ax.plot(positions, wave, label=f"Binary: {binary}", color="blue")
    ax.set_title(f"Binary: {binary}")
    ax.set_xlabel("Cysteine Position")
    ax.set_ylabel("Amplitude")
    ax.set_xticks([1, 2, 3])
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="upper right")

# Adjust layout and show the plot
plt.tight_layout()
plt.show()
