"""
Plots particular solutions to the food chain system
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# ODE System
def lotka(t, vector, a, b, g, d, m, n, p):
    """
    The system to solve, with parameters:
    alpha(a), beta(b), gamma(g), delta(d),
    mu(m), n, and p    
    """
    x,y,z = vector    # unpack input vector
    return a*x-b*x*y, d*x*y-g*y-m*z*y, n*y*z-p*z


def jac(t, vector, a, b, g, d, m, n, p):
    """
    Jacobian
    """
    x,y,z = vector    # unpack input vector
    return np.array([
        [a-b*y, -b*x, 0],
        [d*y, d*x - g - m*z, -m*y],
        [0, n*z, n*y - p]
    ])


# Solution parameters
t0 = 0
tf = 20
pts = 8000

# Set up graph
fig = plt.figure(figsize=[10,14])
ax1 = fig.add_subplot(3, 1, 1, projection='3d')
ax2 = fig.add_subplot(3, 1, 2, projection='3d')
ax3 = fig.add_subplot(3, 1, 3, projection='3d')



# Chosen parameters for system (in order: a, b, g, d, m, n, p)
p1s = [0.5,1,1,1,1,1,1]     # ab < pn
p2s = [1,1,1,1,1,1,1]       # ab = pn
p3s = [1.1,1,1,1,1,1,1]     # ab > pn

i = 0
for pair in [(ax1,p1s), (ax2, p2s), (ax3, p3s)]:
    # unpack parameter values
    ax, ps = pair
    a, b, g, d, m, n, p = ps[0], ps[1], ps[2], ps[3], ps[4], ps[5], ps[6]

    # Solve system for 2 different ICs
    for IC in [(0.1, 0.1, 0.1), (2,2,2), (3,1,3)]:
        # Unpack ICs
        x0, y0, z0 = IC

        # Solve system and plot solution
        soln = solve_ivp(lotka, (t0,tf), (x0,y0,z0), method="LSODA", t_eval=np.linspace(t0,tf,pts), args=(a,b,d,g,m,n,p), jac=jac)
        xs = soln.y[0]
        ys = soln.y[1]
        zs = soln.y[2]
        ax.plot(xs, ys, zs, label=f"({x0}, {y0}, {z0})")

        # plot fixed points
        ax.scatter([0], [0], [0], c="red")
        if(ps in [p1s, p3s]):
            ax.scatter([g/d], [a/b], [0], c="red")
        else:
            non_isolated_xs = np.linspace(0, 6, 400)
            ax.plot(non_isolated_xs, a/b*np.ones_like(non_isolated_xs), 1/m*(d*non_isolated_xs - g), "red")

    # Logic for setting title of figure
    if (a*n < b*p):
        title = r"$\alpha n < \beta p$"
    elif(a*n == b*p):
        title = r"$\alpha n = \beta p$"
    else:
        title = r"$\alpha n > \beta p$"
    
    ax.set_title(title + f" --- Parameters: {ps}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend(loc="upper left")
    i += 1

fig.tight_layout(h_pad=4)
plt.show()