from numpy import exp, sqrt, pi
from scipy.integrate import nquad, quad
from scipy.constants import hbar, m_e, epsilon_0, e

# Bohr radius:
a = 5.29E-11
# Distance between atoms:
L = .106E-9

# Upper and lower bounds (around 100 times the distance between the atoms):
lower = -1E-8
upper = 1E-8

# Wavefunctions (normalized):
def psi_bonding(x, y, z):
    return .562*( 1/sqrt(pi)*a**(-3/2)*exp(-sqrt(x**2 + y**2 + z**2)/a) + 1/sqrt(pi)*a**(-3/2)*exp(-sqrt((x - L)**2 + y**2 + z**2)/a) )

def psi_antibonding(x, y, z):
    return 1.10*( 1/sqrt(pi)*a**(-3/2)*exp(-sqrt(x**2 + y**2 + z**2)/a) - 1/sqrt(pi)*a**(-3/2)*exp(-sqrt((x - L)**2 + y**2 + z**2)/a) )

# Probability density functions (normalized):
def P_bonding(x, y, z):
    return psi_bonding(x, y, z)**2

def P_antibonding(x, y, z):
    return psi_antibonding(x, y, z)**2

# The electron's potential energy function:
def U(x, y, z):
    # Prevent division by zero (we are basically ignoring the singularities directly on top of the protons):
    if (x == 0 and y == 0 and z == 0) or (x == L and y == 0 and z == 0):
        return 0
    else:
        return -e**2/(4*pi*epsilon_0)*( 1/sqrt(x**2 + y**2 + z**2) + 1/sqrt((x - L)**2 + y**2 + z**2) )

# Laplacians of the wavefunctions (retrieved from laplacian.py):
def laplacian_psi_bonding(x, y, z):
    if (x == 0 and y == 0 and z == 0) or (x == L and y == 0 and z == 0):
        return 0
    else:
        return 0.562*(x**2*exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**3.5*(x**2 + y**2 + z**2)) + (L - x)**2*exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**3.5*(y**2 + z**2 + (L - x)**2)) + x**2*exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**2.5*(x**2 + y**2 + z**2)**(3/2)) + (L - x)**2*exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**2.5*(y**2 + z**2 + (L - x)**2)**(3/2)) - exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**2.5*sqrt(y**2 + z**2 + (L - x)**2)) - exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**2.5*sqrt(x**2 + y**2 + z**2)))/sqrt(pi) + 0.562*(y**2*exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**3.5*(y**2 + z**2 + (L - x)**2)) + y**2*exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**3.5*(x**2 + y**2 + z**2)) + y**2*exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**2.5*(y**2 + z**2 + (L - x)**2)**(3/2)) + y**2*exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**2.5*(x**2 + y**2 + z**2)**(3/2)) - exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**2.5*sqrt(y**2 + z**2 + (L - x)**2)) - exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**2.5*sqrt(x**2 + y**2 + z**2)))/sqrt(pi) + 0.562*(z**2*exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**3.5*(y**2 + z**2 + (L - x)**2)) + z**2*exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**3.5*(x**2 + y**2 + z**2)) + z**2*exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**2.5*(y**2 + z**2 + (L - x)**2)**(3/2)) + z**2*exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**2.5*(x**2 + y**2 + z**2)**(3/2)) - exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**2.5*sqrt(y**2 + z**2 + (L - x)**2)) - exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**2.5*sqrt(x**2 + y**2 + z**2)))/sqrt(pi)

def laplacian_psi_antibonding(x, y, z):
    if (x == 0 and y == 0 and z == 0) or (x == L and y == 0 and z == 0):
        return 0
    else:
        return 1.1*(x**2*exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**3.5*(x**2 + y**2 + z**2)) - (L - x)**2*exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**3.5*(y**2 + z**2 + (L - x)**2)) + x**2*exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**2.5*(x**2 + y**2 + z**2)**(3/2)) - (L - x)**2*exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**2.5*(y**2 + z**2 + (L - x)**2)**(3/2)) + exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**2.5*sqrt(y**2 + z**2 + (L - x)**2)) - exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**2.5*sqrt(x**2 + y**2 + z**2)))/sqrt(pi) + 1.1*(-y**2*exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**3.5*(y**2 + z**2 + (L - x)**2)) + y**2*exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**3.5*(x**2 + y**2 + z**2)) - y**2*exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**2.5*(y**2 + z**2 + (L - x)**2)**(3/2)) + y**2*exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**2.5*(x**2 + y**2 + z**2)**(3/2)) + exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**2.5*sqrt(y**2 + z**2 + (L - x)**2)) - exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**2.5*sqrt(x**2 + y**2 + z**2)))/sqrt(pi) + 1.1*(-z**2*exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**3.5*(y**2 + z**2 + (L - x)**2)) + z**2*exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**3.5*(x**2 + y**2 + z**2)) - z**2*exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**2.5*(y**2 + z**2 + (L - x)**2)**(3/2)) + z**2*exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**2.5*(x**2 + y**2 + z**2)**(3/2)) + exp(-sqrt(y**2 + z**2 + (L - x)**2)/a)/(a**2.5*sqrt(y**2 + z**2 + (L - x)**2)) - exp(-sqrt(x**2 + y**2 + z**2)/a)/(a**2.5*sqrt(x**2 + y**2 + z**2)))/sqrt(pi)

# Compute the energy expectation value (the integral over psi*-h_bar**2/2m*laplacian(psi) + psi*U*psi):
def E_bonding(x, y, z):
    return psi_bonding(x, y, z)*( -hbar**2/(2*m_e)*( laplacian_psi_bonding(x, y, z) ) + U(x, y, z)*psi_bonding(x, y, z) )

def E_antibonding(x, y, z):
    return psi_antibonding(x, y, z)*( -hbar**2/(2*m_e)*( laplacian_psi_antibonding(x, y, z) ) + U(x, y, z)*psi_antibonding(x, y, z) )

# Compute the integrals using 1000 steps and then convert the energy expectation values in Joules to eV:
integral_bonding = nquad(E_bonding, [[lower, upper], [lower, upper], [lower, upper]])

print("Energy Expectation Value of the Bonding State (eV): ")
print(integral_bonding[0]*6.24E18)

print()

integral_antibonding = nquad(E_antibonding, [[lower, upper], [lower, upper], [lower, upper]])

print("Energy Expectation Value of the Antibonding State (eV): ")
print(integral_antibonding[0]*6.24E18)