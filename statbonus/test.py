import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from collections import defaultdict

# Directions for square lattice
directions = {
    'R': np.array([1, 0]),
    'L': np.array([-1, 0]),
    'U': np.array([0, 1]),
    'D': np.array([0, -1])
}

# All 3-step walks (4^3 = 64)
all_walks = list(product(directions.keys(), repeat=3))

# Count occurrences of endpoints
regular_counts = defaultdict(int)
self_avoiding_counts = defaultdict(int)

for walk in all_walks:
    pos = np.array([0, 0])
    visited = {tuple(pos)}
    self_avoiding = True

    for step in walk:
        pos = pos + directions[step]
        if tuple(pos) in visited:
            self_avoiding = False
        visited.add(tuple(pos))

    endpoint = tuple(pos)
    regular_counts[endpoint] += 1
    if self_avoiding:
        self_avoiding_counts[endpoint] += 1

# Convert to arrays for plotting
def plot_endpoint_counts(counts, title):
    # Determine bounds
    all_points = list(counts.keys())
    xs, ys = zip(*all_points)
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)

    # Prepare grid
    grid_width = x_max - x_min + 1
    grid_height = y_max - y_min + 1
    grid = np.zeros((grid_height, grid_width))

    for (x, y), count in counts.items():
        grid[y_max - y, x - x_min] = count  # flip y for plotting

    # Plot
    plt.figure(figsize=(6, 6))
    plt.imshow(grid, cmap='viridis', interpolation='nearest')
    plt.colorbar(label='Number of Walks Ending at Point')
    plt.title(title)
    plt.xticks(ticks=np.arange(grid_width), labels=np.arange(x_min, x_max + 1))
    plt.yticks(ticks=np.arange(grid_height), labels=range(y_max, y_min - 1, -1))
    plt.xlabel("x")
    plt.ylabel("y")

    # Annotate counts
    for (x, y), count in counts.items():
        plt.text(x - x_min, y_max - y, str(count),
                 ha='center', va='center', color='white', fontsize=8)

    plt.grid(False)
    plt.tight_layout()
    plt.show()

# Plot regular and self-avoiding endpoint distributions
plot_endpoint_counts(regular_counts, "Regular 3-Step Random Walk Endpoint Counts")
plot_endpoint_counts(self_avoiding_counts, "Self-Avoiding 3-Step Random Walk Endpoint Counts")

