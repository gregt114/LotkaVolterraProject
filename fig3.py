"""
Plots the phase plane to the improved
Lotka-Volterra system
"""


# Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Create figure and axes
fig = plt.figure(figsize=[8,6])
ax  = fig.add_subplot(1,1,1)

# ODE System
def lotka(t, vector, a, b, g, d, k):
    """
    The standard Lotka-Volterra system with
    parameters alpha(a), beta(b), gamma(g), and delta(d).
    """
    x,y = vector    # unpack input vector
    return a*x*(1-x/k) - b*x*y, d*x*y - g*y


# Chosen parameters for system
a = 1
b = 0.1
g = 9
d = 1
k = 4


# Phase Plane parameters
xbound = 15
ybound = 15
n      = 800


# Solve system
xs,ys = np.meshgrid(np.linspace(0, xbound, n), np.linspace(0, ybound, n))
dX,dY = lotka(0, (xs,ys), a,b,g,d,k)  # X,Y derivatives at each point

# Plotting
ax.streamplot(xs, ys, dX, dY, density=2.5, color=np.sqrt(dX**2 + dY**2), cmap="inferno")
ax.scatter([0, k, g/d], [0, 0, a/b*(1-g/(d*k))], c="red", s=34)

# Plotting
ax.set_title(r"ILV System for $\mathbf{\delta k - \gamma < -4 \alpha}$  " +   f"(a={a}, b={b:.1f}, g={g}, d={d}, k={k})")
ax.set_xlabel("Prey")
ax.set_ylabel("Predators")
ax.set_xlim(-0.5, xbound)
ax.set_ylim(-0.5, ybound)

ax.grid()

# plt.show()

# Save plot to image file
plt.savefig("images/fig3.png")
