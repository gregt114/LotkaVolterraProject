"""
Plots a particular solution to the standard
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

# Inital conditions
x0 = 0.9
y0 = 0.9

# Solution domain
t0 = 0
tf = 30
n  = 3000   # number of points in integration

# Solve system
sol = solve_ivp(lotka, (t0,tf), (x0,y0), method="RK45", t_eval=np.linspace(t0,tf,n), args=(a,b,g,d))

# Plotting
ax.plot(sol.t, sol.y[0], label="Prey")
ax.plot(sol.t, sol.y[1], label="Predators")

ax.set_title(f"Lotka Voltera System: a={a:.2f}, b={b:.2f}, g={g}, d={d}")
ax.set_xlabel("Time")
ax.set_ylabel("Count")


ax.legend()
ax.grid()

# plt.show()

# Save plot to image file
plt.savefig("images/fig1.png")
