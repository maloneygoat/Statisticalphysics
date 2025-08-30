import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

def P_reduced(Vr, Tr):
    return (8 * Tr) / (3 * Vr - 1) - 3 / Vr**2

def find_coexistence_volumes_and_pressure(Tr):
    Vr_vals = np.linspace(0.35, 3.0, 1000)
    Pr_vals = P_reduced(Vr_vals, Tr)

    Pr_min = np.nanmin(Pr_vals)
    Pr_max = np.nanmax(Pr_vals)

    def volume_diff(P_guess):
        roots = []
        for i in range(len(Vr_vals) - 1):
            if (Pr_vals[i] - P_guess) * (Pr_vals[i+1] - P_guess) < 0:
                root = brentq(lambda Vr: P_reduced(Vr, Tr) - P_guess, Vr_vals[i], Vr_vals[i+1])
                roots.append(root)
        if len(roots) == 2:
            V1, V2 = roots
            Vr_range = np.linspace(V1, V2, 500)
            f = P_reduced(Vr_range, Tr) - P_guess
            area = np.trapz(f, Vr_range)
            return area
        return 1e3

    try:
        P_eq = brentq(volume_diff, Pr_min + 0.01, Pr_max - 0.01)
        coexisting_volumes = []
        for i in range(len(Vr_vals) - 1):
            if (Pr_vals[i] - P_eq) * (Pr_vals[i+1] - P_eq) < 0:
                root = brentq(lambda Vr: P_reduced(Vr, Tr) - P_eq, Vr_vals[i], Vr_vals[i+1])
                coexisting_volumes.append(root)
        if len(coexisting_volumes) == 2:
            return coexisting_volumes[0], coexisting_volumes[1], P_eq
    except ValueError:
        pass
    return None, None, None

# Compute binodal points
Tr_list = np.linspace(0.6, 0.99, 40)
Vr_liq_list = []
Vr_vap_list = []
Pr_list = []

for Tr in Tr_list:
    Vr_liq, Vr_vap, Pr_eq = find_coexistence_volumes_and_pressure(Tr)
    if Vr_liq and Vr_vap:
        Vr_liq_list.append(Vr_liq)
        Vr_vap_list.append(Vr_vap)
        Pr_list.append(Pr_eq)

# Add the critical point
Vr_liq_list.append(1.0)
Vr_vap_list.append(1.0)
Pr_list.append(1.0)

# Plot binodal line in P_r vs V_r
plt.figure(figsize=(10, 6))
plt.plot(Vr_liq_list, Pr_list, label="Liquid branch", color="blue")
plt.plot(Vr_vap_list, Pr_list, label="Vapor branch", color="orange")
plt.scatter([1.0], [1.0], color="red", label="Critical Point")

plt.xlabel("Reduced Volume $V_r = V/V_c$")
plt.ylabel("Reduced Pressure $P_r = P/P_c$")
plt.title("Binodal Curve in $P_r$ vs. $V_r$ (Van der Waals Fluid - Mercury)")
plt.legend()
plt.grid(True)
plt.xlim(0.35, 3)
plt.ylim(0, 1.2)
plt.show()