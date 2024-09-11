import matplotlib.pyplot as plt
import numpy as np

# Function to generate Pascal's triangle up to n rows
def generate_pascals_triangle(n):
    triangle = [[1]]
    for i in range(1, n):
        row = [1]
        for j in range(1, i):
            row.append(triangle[i-1][j-1] + triangle[i-1][j])
        row.append(1)
        triangle.append(row)
    return triangle

# Function to plot Pascal's triangle
def plot_pascals_triangle(triangle):
    fig, ax = plt.subplots()
    ax.axis('off')  # Turn off the axis

    # Create a grid of numbers in Pascal's triangle
    max_rows = len(triangle)
    for i, row in enumerate(triangle):
        for j, num in enumerate(row):
            ax.text(j - i / 2, -i, str(num), ha='center', va='center')

    ax.set_aspect('equal')
    plt.xlim(-max_rows // 2, max_rows // 2)
    plt.ylim(-max_rows, 0)

    # Save the image at 300 DPI
    plt.savefig("pascals_triangle.png", dpi=300, bbox_inches='tight')
    plt.show()

# Generate Pascal's triangle with 11 rows
pascals_triangle = generate_pascals_triangle(11)

# Plot and export the triangle
plot_pascals_triangle(pascals_triangle)
