# Install the necessary libraries 

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations


# Initialize parameters
num_cysteines = 10  # Total number of cysteines

# Generate all unique proteoforms
proteoforms = []
for i in range(num_cysteines + 1):  # From 0 to 10 cysteines oxidized
    for comb in combinations(range(num_cysteines), i):
        proteoform = np.zeros(num_cysteines, dtype=int)
        proteoform[list(comb)] = 1
        proteoforms.append(proteoform)

# Convert list of arrays to a 2D NumPy array
data = np.array(proteoforms)

# Create a heatmap
plt.figure(figsize=(12, 60))
ax = sns.heatmap(data, cmap="Greys", cbar=False, linewidths=0.5, linecolor='gray')
ax.set_facecolor('black')  # Set the background color of the cells that have no line
plt.title("Cysteine Redox Proteoforms", fontsize=20)
plt.xlabel("Cysteine Sites", fontsize=15)
plt.ylabel("Proteoforms", fontsize=15)
plt.xticks(ticks=np.arange(0.5, num_cysteines + 0.5), labels=np.arange(1, num_cysteines + 1), fontsize=12)
plt.yticks(fontsize=10)

# Saving the figure with high resolution
plt.savefig('/content/Cysteine_Redox_Proteoforms.png', dpi=300, bbox_inches='tight')
plt.show()


# Provide the link to download the saved file
print("Download your high-resolution heatmap here: [Cysteine Redox Proteoforms](sandbox:/mnt/data/Cysteine_Redox_Proteoforms.png)")
