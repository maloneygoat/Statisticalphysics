import numpy as np
from random import randint
import matplotlib
import matplotlib.pyplot as plt
from collections import deque
np.random.seed(2015)

J = 1.0
k_B = 1.0
Tc = 2 / np.log(1 + np.sqrt(2))  # Critical temperature
temperatures = np.linspace(0.01, 1.2 * Tc, 10000)  # Avoid T=0 to prevent division by zero

def dice_samples(trials):
    prob = {1: 1/2, 2: 1/4, 3: 1/8, 4: 1/16, 5: 1/32, 6: 1/32}
    samples = np.zeros(trials + 1, dtype=int)
    samples[0] = 1
    for i in range(trials):
        a = samples[i]
        b = randint(1, 6) # uniform a priori distribution
        pa = prob[a]
        pb = prob[b]
        if pb >= pa or np.random.rand() < pb / pa:
            samples[i + 1] = b
        else:
            samples[i + 1] = a
    return samples

def summarize(samples):
    '''
    Return the percentage of every face in the samples
    '''
    num_samples = len(samples)
    distribution = {i: (samples == i).sum() * 100 / num_samples for i in [1, 2, 3, 4, 5 ,6]} # 百分数
    return distribution

samples = dice_samples(1000000)
ns = np.array(np.logspace(1, 6, num=50), dtype=int)

distributions = {i: np.zeros(50) for i in [1, 2, 3, 4 ,5 ,6]}
for index in range(50):
    n = ns[index]
    distribution = summarize(samples[:n])
    for i in [1, 2, 3, 4, 5, 6]:
        distributions[i][index] = distribution[i]
for i in [1, 2, 3, 4, 5, 6]:
    plt.plot(ns, distributions[i], label='Face {}'.format(i))
plt.xlabel('MC iterations')
plt.ylabel('Percentage')
plt.ylim(0, 100)
plt.xlim(10, 10 ** 6)
plt.semilogx()
plt.legend()
plt.grid()
plt.title("The Metropolis-Hastings dice")
plt.show()


# First, we need helper function to transform between the (i, j) coordinate of a 2D lattice and a serial one
# a: length of the square lattice's side
def flatten_2d(i, j, a):
    return i * a + j # serial No. = row No. * lenght + colum No.
def unflatten_2d(n, a):
    j = n % a
    i = (n - j) // a
    return i, j

# Generate the adjacency list
def gen_neighbors_1d(N):
    neighbors = np.zeros((N, 2), dtype=int)
    for n in range(N):
        neighbors[n][0] = (n - 1) % N # left
        neighbors[n][1] = (n + 1) % N # right
    return neighbors

def gen_neighbors_2d(a):
    neighbors = np.zeros((a*a, 4), dtype=int)
    for n in range(a*a):
        i, j = unflatten_2d(n, a)
        neighbors[n][0] = flatten_2d(i, (j - 1) % a, a) # left
        neighbors[n][1] = flatten_2d(i, (j + 1) % a, a) # right
        neighbors[n][2] = flatten_2d((i - 1) % a, j, a) # up
        neighbors[n][3] = flatten_2d((i + 1) % a, j, a) # down
    return neighbors


def MH_single_flip(neighbors_list, T, iterations):
    '''
    This function performs single flip MC iterations for an Ising system with arbitrary topology, 
    given by the adjaceny list `neighbors_list`.
    The inital state is chosen randomly.
    
    Returns
    =======
    `magnetization`: magnetization (average molecular spin) at each MC step
    `energy`: total energy of the system at each MC step
    '''
    # Initialization
    size = neighbors_list.shape[0]
    spins = np.random.random_integers(0, 1, size)
    spins[spins == 0] = -1
    # Allocation
    magnetization = np.zeros(iterations + 1)
    energy = np.zeros(iterations + 1)
    magnetization[0] = spins.sum()
    energy[0] = -spins.dot(spins[neighbors_list].sum(axis=1)) / 2
    
    for step in range(iterations):
        n = np.random.randint(0, size) # Choose next state according to the a priori distribution
        delta_E = 2 * spins[n] * spins[neighbors_list[n]].sum()
        if delta_E < 0 or np.random.rand() < np.exp(-delta_E / T):
            # Acceptance
            spins[n] = -spins[n]
            magnetization[step + 1] = magnetization[step] + 2 * spins[n]
            energy[step + 1] = energy[step] + delta_E
        else:
            # Rejection
            magnetization[step + 1] = magnetization[step]
            energy[step + 1] = energy[step]
    return magnetization / size, energy


def plot_magnetization(dimension):
    if dimension == 1:
        neighbors_list = gen_neighbors_1d(400)
    elif dimension == 2:
        neighbors_list = gen_neighbors_2d(20)
    else:
        raise ValueError("Unsupported dimension: choose 1 or 2")

    T_list = [0.001 * Tc, 1.0 * Tc]

    for T in T_list:
        magnetization, _ = MH_single_flip(neighbors_list, T, 100000)

        # Create a new figure for each temperature
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle(f"Monte Carlo Simulation (T = {T:.3f}) - {dimension}D Ising Model", fontsize=14)

        # Plot magnetization vs. iteration
        ax1.plot(magnetization, color='blue')
        ax1.set_ylim(-1, 1)
        ax1.set_ylabel('Magnetization')
        ax1.set_xlabel('Iterations')
        ax1.set_title('Time Series')
        ax1.grid(True)

        # Plot histogram of magnetization
        ax2.hist(magnetization, bins=np.linspace(-1, 1, 20), orientation='horizontal', color='orange')
        ax2.set_ylim(-1, 1)
        ax2.set_xlabel('Counts')
        ax2.set_title('Magnetization Distribution')
        ax2.grid(True)

        plt.tight_layout()
        plt.show()

# Example usage
#plot_magnetization(1)
plot_magnetization(2)

############################################


# Cluster MC

def cluster_MC(neighbors_list, T, iterations):
    p = 1 - np.exp(-2 / T) # "magic number"
    # Initialization
    size = neighbors_list.shape[0]
    spins = np.random.random_integers(0, 1, size)
    spins[spins == 0] = -1
    # Allocation
    magnetization = np.zeros(iterations + 1)
    magnetization[0] = spins.sum()
    energy = np.zeros(iterations + 1)
    energy[0] = -spins.dot(spins[neighbors_list].sum(axis=1))
    
    for step in range(iterations):
        # Use a deque to implement breadth-first search
        n0 = np.random.randint(0, size)
        sign = spins[n0]
        cluster = set([n0])
        pockets = deque([n0])
        finished = False
        while not finished:
            try:
                n = pockets.popleft()
                neighbors = neighbors_list[n]
                for neighbor in neighbors:
                    if spins[neighbor] == sign and neighbor not in cluster and np.random.rand() < p:
                        cluster.add(neighbor)
                        pockets.append(neighbor)
            except IndexError:
                finished = True
        # Flip the cluster
        cluster = np.fromiter(cluster, dtype=int)
        spins[cluster] = -sign
        magnetization[step + 1] = magnetization[step] - 2 * sign * len(cluster)
        energy[step + 1] = -spins.dot(spins[neighbors_list].sum(axis=1))
    return magnetization / size, energy / 2 # Every pair is counted two times



def plot_magnetization_cluster(dimension):
    if dimension == 1:
        neighbors_list = gen_neighbors_1d(400)
    elif dimension == 2:
        neighbors_list = gen_neighbors_2d(20)
    else:
        raise ValueError("Unsupported dimension: choose 1 or 2")

    T_list = [0.5, 1.0, 1.5, 1.8, 2.0, 2.2, 2.4, 3.0, 3.5]

    for T in T_list:
        magnetization, _ = cluster_MC(neighbors_list, T, 1000)

        # Create a new figure for each temperature
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
        fig.suptitle(f"Monte Carlo Simulation (T = {T}) - {dimension}D Ising Model", fontsize=14)

        # Time series plot
        ax1.plot(magnetization)
        ax1.set_ylim(-1, 1)
        ax1.set_ylabel('Magnetization')
        ax1.set_xlabel('Iterations')
        ax1.set_title('Time Series')

        # Histogram plot
        ax2.hist(magnetization, bins=np.linspace(-1, 1, num=20), orientation='horizontal')
        ax2.set_ylim(-1, 1)
        ax2.set_xlabel('Counts')
        ax2.set_title('Magnetization Distribution')

        plt.tight_layout()
        plt.show()

# Example usage
plot_magnetization_cluster(2)


#### Actual plot

# The code below runs VERY SLOWLY
points = 100
dimension = 20
iterations = 10000

neighbors_list = gen_neighbors_2d(dimension)
Ts = np.linspace(1.0, 2.1*Tc, points)
Ms = np.zeros(points)
for i in range(points):
    T = Ts[i]
    magnetization, _ = cluster_MC(neighbors_list, T, iterations)
    Ms[i] = np.abs(magnetization).mean()
    print("Iteration for T = {:.3f} complete".format(T)) # Uncomment to print progress as the simulation goes

Onsager_Tc = Tc
Ts_plot = np.linspace(1.0, 2.1*Tc, num=200)
Onsager_Ms = np.zeros(len(Ts_plot))
for i, T in enumerate(Ts_plot):
    if T <= 2.269:
        Onsager_Ms[i] = (1 - (np.sinh(2/T))**(-4))**(1/8)
plt.plot(Ts, Ms, '^', label='simulation')
plt.plot(Ts_plot, Onsager_Ms, '--', label='theoretical')
plt.ylim(0, 1)
plt.legend()
plt.xlabel('Temperature')
plt.ylabel('Magnetization')
plt.grid()
plt.title("Theoretical vs. Numerical values for spontaneous magnitization")
plt.show()