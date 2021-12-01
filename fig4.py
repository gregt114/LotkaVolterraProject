"""
Plots the phase plane to the improved
Lotka-Volterra system
"""


# Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Create figure and axes
fig = plt.figure(figsize=[16,8])
ax1 = fig.add_subplot(1,3,1)
ax2 = fig.add_subplot(1,3,2)
ax3 = fig.add_subplot(1,3,3)

# ODE System
def lotka(t, vector, a, b, g, d, k):
    """
    The Improved Lotka-Volterra system with
    parameters alpha(a), beta(b), gamma(g), and delta(d).
    """
    x,y = vector    # unpack input vector
    return a*x*(1-x/k) - b*x*y, d*x*y - g*y


# Phase Plane parameters
xbound = 8
ybound = 8
n      = 800


# Chosen parameters for system (in order: a, b, g, d, k)
p1s = [1, 1, 1, 1, 2]
p2s = [8, 1, 1, 1, 2]
p3s = [12, 1, 1, 1, 2]

for pair in [(ax1, p1s), (ax2, p2s), (ax3, p3s)]:
    # Unpack data
    ax, ps    = pair
    a,b,g,d,k = ps
    
    # Logic for displaying correct title based on parameter values
    if(a < 4*d*k*(d*k-g)/g):
        title = r"$\mathbf{a < \frac{4\delta k (\delta k -\gamma)}{\gamma}}$" + f"  (a={a}, b={b}, g={g}, d={d}, k={k})"
    elif(a == 4*d*k*(d*k-g)/g):
        title = r"$\mathbf{a = \frac{4\delta k (\delta k -\gamma)}{\gamma}}$" + f"  (a={a}, b={b}, g={g}, d={d}, k={k})"
    else:
        title = r"$\mathbf{a > \frac{4\delta k (\delta k -\gamma)}{\gamma}}$" + f"  (a={a}, b={b:}, g={g}, d={d}, k={k})"


    # Solve system
    xs,ys = np.meshgrid(np.linspace(0, xbound, n), np.linspace(0, ybound, n))
    dX,dY = lotka(0, (xs,ys), a,b,g,d,k)  # X,Y derivatives at each point

    # Plotting
    ax.streamplot(xs, ys, dX, dY, density=2.5, color=np.sqrt(dX**2 + dY**2), cmap="plasma")
    ax.scatter([0, k, g/d], [0, 0, a/b*(1-g/(d*k))], c="red", s=34)     # Fixed points

    ax.set_title(title)
    ax.set_xlabel("Prey")
    ax.set_ylabel("Predators")
    ax.set_xlim(-0.2, xbound)
    ax.set_ylim(-0.2, ybound)

    ax.grid()

# plt.show()

plt.tight_layout()

# Save plot to image file
plt.savefig("images/fig4.png")
