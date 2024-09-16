import math
import pandas as pd

# Target large number M
M = 2.90e16

# Number of proteoforms (keys in the proteoforms dictionary)
proteoforms = {
    'p0': 0,
    'p10': 10,
    'p20': 20,
    'p30': 30,
    'p40': 40,
    'p50': 50,
    'p60': 60,
    'p70': 70,
    'p80': 80,
    'p90': 90,
    'p100': 100
}
num_proteoforms = len(proteoforms)

# Function to estimate the solution space size using combinatorics
def estimate_solution_space_size(num_molecules, num_proteoforms):
    return math.comb(num_molecules + num_proteoforms - 1, num_proteoforms - 1)

# Now, let's find the number of molecules where the solution space size is close to M
def find_molecule_count_for_M(M, num_proteoforms):
    num_molecules = 1
    results = []  # List to store the results
    while True:
        solution_space_size = estimate_solution_space_size(num_molecules, num_proteoforms)
        print(f"With {num_molecules} molecules: Solution space size = {solution_space_size}")
        
        # Append the current molecule count and solution space size to the results list
        results.append({"Molecule Count": num_molecules, "Solution Space Size": solution_space_size})
        
        if solution_space_size >= M:
            return num_molecules, solution_space_size, results
        
        num_molecules += 1

# Run the function to find the number of molecules for which the solution space size is close to M
molecule_count, solution_space_size, results = find_molecule_count_for_M(M, num_proteoforms)

print(f"\nNumber of molecules required for solution space size to reach or exceed M: {molecule_count}")
print(f"Calculated solution space size: {solution_space_size}")

# Save the results to a CSV file
df_results = pd.DataFrame(results)
df_results.to_csv('solution_space_size_results.csv', index=False)

print("Results saved to 'solution_space_size_results.csv'")
