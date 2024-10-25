from numpy import exp, sqrt, pi
from scipy.integrate import nquad, quad

# Bohr radius:
a = 5.29E-11
# Distance between atoms:
L = .106E-9

# Upper and lower bounds (around 100 times the distance between the atoms):
lower = -1E-8
upper = 1E-8

# Probability density functions (without normalization A^2):
def P_bonding(x, y, z):
    return ( 1/sqrt(pi)*a**(-3/2)*exp(-sqrt(x**2 + y**2 + z**2)/a) + 1/sqrt(pi)*a**(-3/2)*exp(-sqrt((x - L)**2 + y**2 + z**2)/a) )**2

def P_antibonding(x, y, z):
    return ( 1/sqrt(pi)*a**(-3/2)*exp(-sqrt(x**2 + y**2 + z**2)/a) - 1/sqrt(pi)*a**(-3/2)*exp(-sqrt((x - L)**2 + y**2 + z**2)/a) )**2

# Compute the integrals:
integral_bonding = nquad(P_bonding, [[lower, upper], [lower, upper], [lower, upper]])
print("Integral of the Bonding State's Probability Density Function (without normalization): ")
print(integral_bonding)

integral_antibonding = nquad(P_antibonding, [[lower, upper], [lower, upper], [lower, upper]])
print("Integral of the Antibonding State's Probability Density Function (without normalization): ") 
print(integral_antibonding)
