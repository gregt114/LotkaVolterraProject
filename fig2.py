"""
Plots phase plane to the standard
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
def lotka(t, vector, a, b, g, d):
    """
    The standard Lotka-Volterra system with
    parameters alpha(a), beta(b), gamma(g), and delta(d).
    """
    x,y = vector    # unpack input vector
    return a*x - b*x*y, d*x*y - g*y


# Chosen parameters for system
a = 2/3
b = 4/3
g = 1
d = 1

# Phase Plane parameters
xbound = 3
ybound = 3
n      = 800


# Solve system
xs,ys = np.meshgrid(np.linspace(0,xbound, n), np.linspace(0, ybound, n))
dX,dY = lotka(0, (xs,ys), a,b,g,d)  # X,Y derivatives at each point

# Plotting
ax.streamplot(xs, ys, dX, dY, density=2, color=np.sqrt(dX**2 + dY**2), cmap="inferno")
ax.scatter([0, g/d], [0, a/b], c="red", s=34)

ax.set_title(f"Lotka Voltera System: a={a:.2f}, b={b:.2f}, g={g}, d={d}")
ax.set_xlabel("Prey")
ax.set_ylabel("Predators")

ax.grid()

# plt.show()

# Save plot to image file
plt.savefig("fig2.png")
