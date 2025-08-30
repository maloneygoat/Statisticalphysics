import numpy as np
import matplotlib.pyplot as plt
from itertools import product

# Step vectors: right, left, up, down
steps = [(1, 0), (-1, 0), (0, 1), (0, -1)]
#steps = [(5, 0), (-5, 0), (0, 5), (0, -5)]
N = 3  # number of steps

def end_position(path):
    x, y = 0, 0
    for dx, dy in path:
        x += dx
        y += dy
    return (x, y)

def is_self_avoiding(path):
    x, y = 0, 0
    visited = {(x, y)}
    for dx, dy in path:
        x += dx
        y += dy
        if (x, y) in visited:
            return False
        visited.add((x, y))
    return True

# Generate all 3-step walks
all_walks = list(product(steps, repeat=N))

# Containers
normal_endpoints = []
normal_R2 = []
saw_endpoints = []
saw_R2 = []

# Analyze each walk
for walk in all_walks:
    x, y = end_position(walk)
    R2 = x**2 + y**2
    normal_endpoints.append((x, y))
    normal_R2.append(R2)
    if is_self_avoiding(walk):
        saw_endpoints.append((x, y))
        saw_R2.append(R2)

# Statistics
Ω_normal = len(normal_R2)
Ω_SAW = len(saw_R2)
mean_R2_normal = np.mean(normal_R2)
mean_R2_SAW = np.mean(saw_R2)

print(f"(b) Number of configurations:")
print(f"  Ω_normal = {Ω_normal}")
print(f"  Ω_SAW = {Ω_SAW}\n")

print(f"(d) Mean squared end-to-end distance:")
print(f"  ⟨R²⟩_normal = {mean_R2_normal:.3f}")
print(f"  ⟨R²⟩_SAW = {mean_R2_SAW:.3f}")

# (c) Plotting end points on the lattice
plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)
x, y = zip(*normal_endpoints)
plt.scatter(x, y, alpha=0.5)
plt.title("Normal Walk Endpoints")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
plt.axis('equal')

plt.subplot(1, 2, 2)
x_saw, y_saw = zip(*saw_endpoints)
plt.scatter(x_saw, y_saw, color='red', alpha=0.6)
plt.title("Self-Avoiding Walk Endpoints")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
plt.axis('equal')

plt.suptitle("End Points of 3-Step Walks")
plt.tight_layout()
plt.show()