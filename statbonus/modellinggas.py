import numpy as np
from random import randint
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import minimize_scalar
from scipy.optimize import curve_fit
from scipy.interpolate import make_smoothing_spline

import numpy as np
import matplotlib.pyplot as plt

n = 1 # number of moles
R = 0.08314 
# Mercury: 

a = 8.2
b = 0.01696 # Bar
Tc = 1750.1


def pressure(V, T, R, a0, b0):
    A = R * T/(V - b0)
    B = a0 /(V** 2)
    
    return A - B



    

Ts = [0.7 * Tc, 0.85 * Tc, 1.0 * Tc, 1.15 * Tc, 1.3 * Tc]
volume = np.linspace(b + 1e-4, 0.3, 1000)
for i in range(len(Ts)):
    plt.plot(volume, pressure(volume, Ts[i], R, a, b), label='{:.2f} K'.format(Ts[i]))

#plt.plot(volume, dp(volume, Tc, R, a, b), label = "Derivative")
#plt.plot(volume, d2p(volume, Tc, R, a, b))


plt.legend()
plt.xlabel('Volume')
plt.ylabel('Pressure')
#plt.ylim(10, 10 ** 6)
plt.yscale('log')
plt.grid()
plt.title("Isotherms of Mercry at different T")
plt.show()



# Constants for mercury
R = 0.08314       # L·bar/mol·K
a = 8.200         # L²·bar/mol²
b = 0.01696       # L/mol

# Critical constants
Vc = 3 * b
Tc = (8 * a) / (27 * R * b)
Pc = a/(27 * b * b)


print(f"Mercury Critical Constants:\nVc = {Vc:.5f} L/mol\nTc = {Tc:.2f} K\nPc = {Pc:.2f} bar")

# Reduced Van der Waals equation
def P_reduced(Vr, Tr):
    return (8 * Tr) / (3 * Vr - 1) - 3 / Vr**2

#def spinodal(Vr):
    #return -24/((3 * Vr - 1) ** 2) + 6/(Vr ** 3)

def spinodal(Vr):
    return 3/(Vr ** 2) - 2/(Vr ** 3)


# Reduced volume range (must be > 1/3 to avoid singularity)
Vr = np.linspace(0.35, 10, 1000)
Trs = [0.7, 0.85, 1.0, 1.15, 1.3]  # Reduced temperatures

# Plot
plt.figure(figsize=(10, 6))
for Tr in Trs:
    Pr = P_reduced(Vr, Tr)
    plt.plot(Vr, Pr, label=f"Tr = {Tr:.2f}")

plt.xlabel("Reduced Volume $V_r = V/V_c$")
plt.ylabel("Reduced Pressure $P_r = P/P_c$")
plt.title("Reduced Van der Waals Isotherms for Mercury")
plt.grid(True)
plt.legend()
plt.ylim(-2, 3)
plt.xlim(0, 3)
#plt.xscale("log")
plt.show()

### Finding Binodal 
from scipy.optimize import brentq, minimize_scalar
import numpy as np

def find_three_intersections(func, line_val, x_range, num_points=2000, tol=1e-5):
    """
    Find 3 distinct intersection points between a curve and a horizontal line.
    """
    x_vals = np.linspace(x_range[0], x_range[1], num_points)
    y_vals = func(x_vals) - line_val

    roots = []
    for i in range(len(x_vals) - 1):
        if y_vals[i] * y_vals[i+1] < 0:
            try:
                root = brentq(lambda x: func(x) - line_val, x_vals[i], x_vals[i+1])
                # Avoid duplicates (numerical noise)
                if all(abs(root - r) > tol for r in roots):
                    roots.append(root)
            except ValueError:
                continue

    if len(roots) == 3:
        return sorted(roots)
    else:
        return 0

def integrate_between(func, t, a, b, num_points=500):
    """
    Numerically integrates a function func(x) between limits a and b.
    """
    x = np.linspace(a, b, num_points)
    y = func(x, t)
    return np.trapz(y, x)

def leftInt(func, temp, list, line):
    I = integrate_between(func, temp, list[0], list[1])
    area = (list[1] - list[0])*line - I
    return area

def righInt(func, temp, list, line):
    I = integrate_between(func, temp, list[1], list[2])
    area = I - (list[2] - list[1])*line
    return area

#Trs = [0.7, 0.85, 1.0, 1.15, 1.3]
Trs = [0.65, 0.7, 0.85, 0.9, 0.99, 1.0, 1.15, 1.3]
def equalareafinder(guesses, xranges):
    points = []
    for i in range(len(Trs)):
        if Trs[i] < 1:
            Ri = 5
            Li = 0
            last  = 1
            func = lambda Vr: P_reduced(Vr, Trs[i])
            max_res = minimize_scalar(lambda Vr: -func(Vr), bounds=(0.4, 20), method='bounded')
            min_res = minimize_scalar(func, bounds=(0.4, 20), method='bounded')
            max_val = -max_res.fun
            min_val = min_res.fun
            iteration = 0
            
            # We are trying to find a guess where the left and right areas are equal
            def objective_function(guess):
                intersects = find_three_intersections(func, guess, xranges[i])
                if intersects != 0:
                    Li = leftInt(P_reduced, Trs[i], intersects, guess)
                    Ri = righInt(P_reduced, Trs[i], intersects, guess)
                    return abs(Ri - Li)  # Minimize the difference between areas
                else:
                    return 1e5  # Return a large value if no intersections

            # Use minimize_scalar to find the best guess
            result = minimize_scalar(objective_function, bounds=(min_val, max_val), method='bounded')

            # Get the final guess and intersection points
            guess = result.x
            intersects = find_three_intersections(func, guess, xranges[i])
            
            if intersects != 0:
                points.append((guess, intersects))
            else:
                print(f"No valid intersections found at Tr = {Trs[i]}, guess = {guess}")
        
    return points

# Example usage
Tr = 0.85
func = lambda Vr: P_reduced(Vr, Tr)
xranges = [(0.4, 6),(0.45, 5), (0.35, 3), (0.4, 3), (0.4, 3)]
testline = 0.53
guesses = [0.4, 0.47, 0.53, 0.6, 0.8]
intersections = find_three_intersections(func, testline, (0.35, 3))

if intersections != 0:
    print("Three intersection points found at Vr =", intersections)
else:
    print("Line does not intersect the curve three times.")

value = equalareafinder(guesses, xranges)
print("Binodal points:", value)


def model_function(x, a, b, c, d):
    """
    Example model function: A quadratic function.
    Modify this based on the type of curve you're fitting.
    """
    return a * x**3 + b * x**2 + c * x + d

def curve_fitting(x_data, y_data, model_func, initial_guess):
    """
    Fits a curve to the data using the provided model function.

    Parameters:
    - x_data: List or array of x-values
    - y_data: List or array of y-values
    - model_func: The model function to fit (e.g., a polynomial function)
    - initial_guess: Initial guess for the parameters of the model function

    Returns:
    - popt: Optimized parameters for the model
    - pcov: Covariance matrix of the optimized parameters
    """
    # Use curve_fit to fit the model function to the data
    popt, pcov = curve_fit(model_func, x_data, y_data, p0=initial_guess)
    
    return popt, pcov


def plot_binodal(binodal_points):
    """
    Plot the binodal points, where the first part is the pressure (Pr) and
    the second part is the volume at that pressure.
    """
    # Unzip the binodal points into Pr (pressure) and Vr (volume)
    pressures, volumes = zip(*binodal_points)
    pres = []
    vol = []
    for i in range(len(pressures)):
        for j in range(3):
            pres.append(pressures[i])
            vol.append(volumes[i][j])
    print(pres)
    print(vol)
    
    # Example usage
    x_data = np.array([0, vol[0], vol[3], vol[6], vol[9], vol[12], 1, vol[-1], vol[11], vol[8], vol[5], vol[2], 10])  # Example x-values
    y_data = np.array([0, pres[0], pres[3], pres[6], pres[9], pres[12], 1, pres[12], pres[9], pres[6], pres[3], pres[0], 0])  # Example y-values, which fit a quadratic curve
    
    # Initial guess for the parameters [a, b, c] of the model function
    initial_guess = [0.1, 0.1, 5, 2]

    # Fit the data using the model function
    popt, pcov = curve_fitting(x_data, y_data, model_function, initial_guess)

    # Extract the fitted parameters
    a, b, c, d = popt

    # Print the fitted parameters
    print(f"Fitted parameters: a = {a}, b = {b}, c = {c}, d = {d}")

            
    spl = make_smoothing_spline(x_data, y_data, lam=0.000000000001)
    # Plotting
    plt.figure(figsize=(10, 6))

    # Plot the isotherms for different Tr values
    for Tr in Trs:
        Vr = np.linspace(0.35, 10, 500)  # Sample a range of Vr values for plotting the isotherms
        Pr = P_reduced(Vr, Tr)  # Calculate the corresponding Pr using P_reduced function
        plt.plot(Vr, Pr, label=f"Tr = {Tr:.2f}")  # Plot the isotherm for this Tr value

    # Add scatter plot for the binodal points
    plt.scatter(vol, pres, color='blue', marker='o', label="Intersection points")
    plt.plot(x_data, spl(x_data), 'b--', label='Binodal curve')
    # Fitted curve

    # Adding labels and title
    plt.xlabel("Reduced Volume $V_r = V/V_c$")
    plt.ylabel("Reduced Pressure $P_r = P/P_c$")
    plt.title("Reduced Van der Waals Isotherms and Binodal Points")
    plt.grid(True)
    plt.legend()

    # Adjust plot limits to better fit the data
    plt.ylim(-2, 3)
    plt.xlim(0, 7)

    # Show the plot
    plt.show()

# Assuming you already have the `value` variable from the previous `equalareafinder`
# If the `value` contains the binodal points, call the plot function
plot_binodal(value)


# Finally plot your og graphs, spinodal and binodals.
Trs = [0.7, 0.85, 1.0, 1.15, 1.3]

def plot_final(binodal_points):
    """
    Plot the binodal points, where the first part is the pressure (Pr) and
    the second part is the volume at that pressure.
    """
    # Unzip the binodal points into Pr (pressure) and Vr (volume)
    pressures, volumes = zip(*binodal_points)
    pres = []
    vol = []
    for i in range(len(pressures)):
        for j in range(3):
            pres.append(pressures[i])
            vol.append(volumes[i][j])
    print(pres)
    print(vol)
    
    # Example usage
    x_data = np.array([0, vol[0], vol[3], vol[6], vol[9], vol[12], 1, vol[-1], vol[11], vol[8], vol[5], vol[2], 8])  # Example x-values
    y_data = np.array([0, pres[0], pres[3], pres[6], pres[9], pres[12], 1, pres[12], pres[9], pres[6], pres[3], pres[0], 0])  # Example y-values, which fit a quadratic curve
         
    spl = make_smoothing_spline(x_data, y_data, lam=0.000000000001)
    # Plotting
    plt.figure(figsize=(10, 6))

    # Plot the isotherms for different Tr values
    
    Vr = np.linspace(0.35, 10, 500)  # Sample a range of Vr values for plotting the isotherms
    for Tr in Trs:
        Vr = np.linspace(0.35, 10, 500)  # Sample a range of Vr values for plotting the isotherms
        Pr = P_reduced(Vr, Tr)  # Calculate the corresponding Pr using P_reduced function
        plt.plot(Vr, Pr, label=f"Tr = {Tr:.2f}")  # Plot the isotherm for this Tr value

    # Add scatter plot for the binodal points

    plt.plot(x_data, spl(x_data), 'b--', label='Binodal')
    
    spin = spinodal(Vr)
    plt.plot(Vr, spin, 'r--', label="Spinodal")
    # Fitted curve

    # Adding labels and title
    plt.xlabel("Reduced Volume $V_r = V/V_c$")
    plt.ylabel("Reduced Pressure $P_r = P/P_c$")
    plt.title("Reduced Van der Waals Isotherms and Binodal Points")
    plt.grid(True)
    plt.legend()

    # Adjust plot limits to better fit the data
    plt.ylim(-2, 3)
    plt.xlim(0, 7)

    # Show the plot
    plt.show()

# Assuming you already have the `value` variable from the previous `equalareafinder`
# If the `value` contains the binodal points, call the plot function
plot_final(value)
