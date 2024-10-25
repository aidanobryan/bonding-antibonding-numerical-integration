from sympy import symbols, exp, sqrt, pi, diff
from sympy.vector import CoordSys3D

# In order to compute the energy expectation value, we must first symbolically compute the Laplacian of the wavefunctions.

# Define the symbols:
x, y, z, a, L = symbols('x y z a L')

# Define the wavefunctions (normalized):
psi_bonding = .562*( 1/sqrt(pi)*a**(-3/2)*exp(-sqrt(x**2 + y**2 + z**2)/a) + 1/sqrt(pi)*a**(-3/2)*exp(-sqrt((x - L)**2 + y**2 + z**2)/a) )

psi_antibonding = 1.10*( 1/sqrt(pi)*a**(-3/2)*exp(-sqrt(x**2 + y**2 + z**2)/a) - 1/sqrt(pi)*a**(-3/2)*exp(-sqrt((x - L)**2 + y**2 + z**2)/a) )

# Compute the Laplacian of the wavefunctions:
def get_laplacian_psi_bonding():
    return diff(psi_bonding, x, x) + diff(psi_bonding, y, y) + diff(psi_bonding, z, z)

print("Laplacian of the Bonding State: ")
print(get_laplacian_psi_bonding())

print()

def get_laplacian_psi_antibonding():
    return diff(psi_antibonding, x, x) + diff(psi_antibonding, y, y) + diff(psi_antibonding, z, z)

print("Laplacian of the Antibonding State: ")
print(get_laplacian_psi_antibonding())
