import itertools
import numpy as np
import pandas as pd

# Proteoforms and their corresponding oxidation states
proteoforms = {
    'alpha': 0,
    'beta': 20,
    'gamma': 40,
    'delta': 60,
    'epsilon': 80,
    'zeta': 100
}

# Target oxidation percentage
target_oxidation = 20

# Fixed number of molecules (e.g., 10 molecules)
num_molecules = 10

# Possible integer counts of molecules in each state (must sum up to num_molecules)
molecule_counts = list(range(num_molecules + 1))

# Function to calculate weighted average
def calculate_weighted_average(proteoform_combination):
    total_molecules = sum(proteoform_combination.values())
    oxidation_sum = sum(proteoforms[p] * (count / total_molecules) for p, count in proteoform_combination.items())
    return oxidation_sum

# Generate all possible combinations of counts summing to num_molecules
combinations = [comb for comb in itertools.product(molecule_counts, repeat=6) if sum(comb) == num_molecules]

# Store valid solutions
valid_solutions = []

# Check each combination for matching the target oxidation level
for combination in combinations:
    proteoform_combination = dict(zip(proteoforms.keys(), combination))
    if np.isclose(calculate_weighted_average(proteoform_combination), target_oxidation):
        valid_solutions.append(proteoform_combination)

# Display the valid solutions
df_solutions = pd.DataFrame(valid_solutions)
print("Possible solutions for 20% oxidation state with fixed number of molecules:")
print(df_solutions)

# Save the solutions to a CSV file (optional)
df_solutions.to_csv('oxidation_state_solutions.csv', index=False)
