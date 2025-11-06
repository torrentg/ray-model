import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import newton
import matplotlib.pyplot as plt

# Constants
G = 1             # Gravitational constant
M = 1             # Mass of the sphere
D = 10            # Distance from the center of the sphere
R = 1             # Radius of the sphere
C = 1             # Speed of light
V = 0.5           # Velocity of the sphere

# Derived constants
s = D / R
beta = V / C
gamma = 1 / math.sqrt(1 - beta ** 2)
rho_0 = M / ((4 / 3) * math.pi * R ** 3)

# Define A(mu), B(mu), and C(mu) for Example 2
def coef_A(mu):
    return 1 + 2 * beta * mu + beta ** 2

def coef_B(mu):
    return D * (mu + beta)

def coef_C(mu):
    return D ** 2 - R ** 2

def discriminant(mu):
    val = coef_B(mu)**2 - coef_A(mu) * coef_C(mu)
    if val < 0:
        return 0
    return val

def L(mu):
    return 2 * R * math.sqrt(discriminant(mu)) / coef_A(mu)

def integrand(mu):
    return L(mu) * (mu + beta) * (1 + beta * mu)

def find_mu0():
    return (-beta - math.sqrt((s * s - 1) * (s * s - beta * beta))) / s ** 2
    #x0 = -1
    #root = newton(discriminant, x0)
    #return root


mu_0 = find_mu0()
print(f"mu_0: {mu_0}")

# Perform the integration
result, _ = quad(integrand, -1, mu_0)

g_ray = -2 * math.pi * G * rho_0 * result / (1 - beta ** 2)
print(f"g_ray: {g_ray}")

g_newton = -G * M / D**2
print(f"g_newton: {g_newton}")

g_rel = g_newton * (1 - beta**2) / (1 + beta) ** 3
print(f"g_rel: {g_rel}")

#---

x_values = []
y_values_ray = []
y_values_rel = []

for i in range(0, 100):

    V = C * i / 100.0
    #D = R + i
    #R = 0.1 + i * 0.2
    x_values.append(V)

    s = D / R
    beta = V / C
    gamma = 1 / math.sqrt(1 - beta ** 2)
    rho_0 = M / ((4 / 3) * math.pi * R ** 3)

    mu_0 = find_mu0()
    result, _ = quad(integrand, -1, mu_0)
    g_ray =  2 * math.pi * G * rho_0 * result / (1 - beta ** 2)
    y_values_ray.append(g_ray)

    g_rel = (-G * M / D**2) * (1 - beta**2) / (1 + beta) ** 3
    y_values_rel.append(g_rel)

# Create the plot
plt.plot(x_values, y_values_ray, label='g_ray', color='blue', linestyle='-') #, marker='o')
plt.plot(x_values, y_values_rel, label='g_rel', color='red', linestyle='-')#, marker='x')

# Add labels and title
plt.xlabel('Velocity (V)')
plt.ylabel('Gravitational Field (g)')
plt.title('Gravitational Field vs Velocity')

# Add a legend
plt.legend()

# Show the plot
plt.grid(True)
plt.show()
