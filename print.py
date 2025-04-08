import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
from matplotlib.colors import SymLogNorm
from scipy.optimize import curve_fit
import re

import scipy.constants as const
from scipy.integrate import simps

# File path
file_path = "./chemistry/swarm.TXT"

with open(file_path, "r") as file:
    lines = file.readlines()
# Locate the "Rate coefficients" section
rate_section = False
transport_section = False
inv_rate_coef = False
E_N_list = []
C2_list = []
C3_list = []
inv_C2 = []

Mobility_list = []
Difusion_list = []

for line in lines:
    if "Rate coefficients (m3/s)" in line:
        rate_section = True

    elif rate_section and re.match(r"^\s*\d+", line):  # Matches data rows
        parts = line.split()
        if len(parts) >= 5:
            E_N_list.append(float(parts[2]))  # Convert to float
            C2_list.append(float(parts[4]))
            C3_list.append(float(parts[5]))
    elif rate_section and line.strip() == "":  # Stop at the next blank line
        break


for line in lines:
    if "Transport coefficients" in line:
        transport_section = True

    elif transport_section and re.match(r"^\s*\d+", line):  # Matches data rows
        parts = line.split()
        if len(parts) >= 6:
            #E_N_list.append(float(parts[1]))  # Convert to float
            Mobility_list.append(float(parts[3]))
            Difusion_list.append(float(parts[4]))
    elif transport_section and line.strip() == "":  # Stop at the next blank line
        break

for line in lines:
    if "Inverse rate coefficients" in line:
        inv_rate_coef = True

    elif inv_rate_coef and re.match(r"^\s*\d+", line):  # Matches data rows
        parts = line.split()
        if len(parts) >= 6:
            #E_N_list.append(float(parts[1]))  # Convert to float
            inv_C2.append(float(parts[4]))
    elif inv_rate_coef and line.strip() == "":  # Stop at the next blank line
        break



# Print lists
print("E/N List:", E_N_list) # Now its electron energy
print("C2 List:", C2_list)
print("C3 List:", C3_list)
print("InvC2 List: ", inv_C2)

print("E/N List:", E_N_list)
print("Mobility List:", Mobility_list)
print("Difusion List:", Difusion_list)

def arctan(x,a,b,c,d): # For the mobility
    return a * np.arctan(np.log10(x) * b - c) + d

def line(x,a,b):
    return a*np.exp(b*np.log(x))

def b(x,a,b,c,d): # For the Ionization Rates

    #return 4.1e-14  * (np.exp(-0.5*( np.log(x*0.00015) / 1.3)**2)) #C2
    #return 1.4e-13  * (np.exp(-0.5*( np.log(x*0.00008) / 1.2)**2)) #C3
    #return a  * (np.exp(-0.5*( np.log(x*b) / c)**2)) # General

    #return a  * (np.exp(b*( np.log(x)/c)**2)) # General
    return a  * (np.exp(-0.5*( (np.log(x*b)) / c)**2)) # General

def cubic(x,a,b,c,d,e):

    return a*np.exp(b*np.log(x)**3 + c*np.log(x)**2+d*np.log(x)) + e
    
# Fit the data with free parameters
poptC2, pcov = curve_fit(
    b, 
    E_N_list[2:], 
    C2_list[2:], 
    p0=[ 1.4e-13, 0.00008, 1.6,0]  # Initial guesses for a, b, and c
)

# Fit the data with free parameters
poptC3, pcov = curve_fit(
    b, 
    E_N_list[2:], 
    C3_list[2:], 
    p0=[1.4e-13, 0.00008, 1.6,0]  # Initial guesses for a, b, and c
)

poptinvC2, pcov = curve_fit(
    b, 
    E_N_list, 
    inv_C2, 
    p0=[1.4e-13, 0.00008, 1.6,0]  # Initial guesses for a, b, and c
)

# Extract fitted parameters
a_fitC2, b_fitC2, c_fitC2,d = poptC2
a_fitC3, b_fitC3, c_fitC3,d3 = poptC3
a_invC2, b_invC2, c_invC2,d_invC2 = poptinvC2

print("C2 optimal values: ", poptC2)
print("C3 optimal values: ", poptC3)
print("Inverse C2 optimal values: ", poptinvC2)
x = np.logspace(-3, 4,100)

plt.plot(x, b(x, a_invC2, b_invC2, c_invC2,d_invC2))
plt.scatter(E_N_list, inv_C2)
plt.plot()
plt.xscale("log")
#plt.yscale("log")

plt.show()

plt.plot(x, b(x, a_fitC2, b_fitC2, c_fitC2,d))
plt.scatter(E_N_list, C2_list)
plt.plot()
plt.xscale("log")
#plt.yscale("log")

plt.show()
#print(b(x, a_fitC3, b_fitC3, c_fitC3,d3))
plt.plot(x, b(x, a_fitC3, b_fitC3, c_fitC3,d3))
plt.scatter(E_N_list, C3_list)
plt.plot()
plt.xscale("log")
#plt.yscale("log")

plt.show()

# Fit the data with free parameters
popt, pcov = curve_fit(
    line, 
    E_N_list[7:17], 
    Mobility_list[7:17], 
    p0=[-1e25, 0]
)
#arctan,
#p0=[-1e25, 1, -5,1e25]  # Initial guesses for a, b, and c

# Extract fitted parameters
#a_fit, b_fit, c_fit, d_fit = popt
a_fit, b_fit = popt
#print(a_fit, b_fit, c_fit, d_fit)
print(a_fit, b_fit)

#print(f)
plt.plot(x,line(x, a_fit,b_fit), #arctan(x, a_fit, b_fit, c_fit, d_fit), 
         label=f'Fitted arctan: a={a_fit:.2e}, b={b_fit:.2e}', 
         color='red')
plt.scatter(E_N_list[7:17], Mobility_list[7:17])
plt.plot()
plt.xscale("log")
plt.yscale("log")
plt.show()

# Fit the data with free parameters
poptcubic, pcov = curve_fit(
    cubic, 
    E_N_list[2:], 
    Difusion_list[2:], 
    p0=[ 1.4e25, 1e-1, -6e-1,-3e-1 , -3e-0],  # Initial guesses for a, b, and c
    maxfev=5000
)

a_fit, b_fit, c_fit, d_fit, e_fit = poptcubic
print("Cubic: ",poptcubic)
plt.plot(x, cubic(x,a_fit, b_fit, c_fit, d_fit, e_fit))
plt.scatter(E_N_list[7:17], Difusion_list[7:17])
plt.plot()
plt.xscale("log")
plt.yscale("log")
plt.show()

# File path
file_path2 = "./chemistry/Cross section.txt"

with open(file_path2, "r") as file2:
    lines2 = file2.readlines()
# Locate the "Rate coefficients" section
ionization = False

electron_energy = []
cross_section = []

for line in lines2:
    print(line)
    if "IONIZATION" in line:
        ionization = True

    elif rate_section and re.match(r"^\s*\d+", line):  # Matches data rows
        parts = line.split()
        if len(parts) >= 2:
            electron_energy.append(float(parts[0]) - 11.5)  # Convert to float
            cross_section.append(float(parts[1]))
    elif rate_section and line.strip() == "":  # Stop at the next blank line
        break

print(electron_energy)
print(cross_section)

# Constants
m_e = const.electron_mass  # Electron mass (kg)
eV_to_J = const.e  # Conversion from eV to Joules

print(eV_to_J)

# Function to compute electron velocity (m/s) from energy (eV)
def electron_velocity(E):
    E = np.array(E, dtype=float)  # Ensure E is a NumPy array
    return np.sqrt(2 * E * eV_to_J / m_e)

# Function to compute Maxwellian EEDF (normalized)
def maxwellian_eedf(E, Te):
    E = np.array(E, dtype=float)  # Ensure E is a NumPy array
    kB_Te = Te*2/3  # Temperature in eV (since kB * Te in eV units)
    norm_factor = 2 / (np.sqrt(np.pi) * (kB_Te)**(3/2))
    return norm_factor * np.sqrt(E) * np.exp(-E / (kB_Te))

# Compute ionization rate for a given electron temperature
def ionization_rate(Te):
    v = electron_velocity(electron_energy)  # Compute electron velocities (array)
    f = maxwellian_eedf(electron_energy, Te)  # Compute Maxwellian EEDF (array)
    integrand = cross_section * v * f  # Compute element-wise product

    # Numerical integration using Simpson's rule
    R = simps(integrand, electron_energy)
    return R  # Ionization rate (m^3/s)

C5 = []
a = 0

for i in E_N_list:
    C5.append(ionization_rate(i))

    print("Calc here for energy: ", C5[a], "Bolsig online - ", C3_list[a])
    a+=1;

print(C5)

poptC5, pcov = curve_fit(
    b, 
    E_N_list, 
    C5, 
    p0=[1.4e-13, 0.00008, 1.6,0]  # Initial guesses for a, b, and c
)

a_C5, b_C5, c_C5,d_C5 = poptC5
print("Step Ionization optimal values: ", poptC5)

plt.plot(x, b(x, a_C5, b_C5, c_C5,d_C5))
plt.scatter(E_N_list, C5)
plt.plot()
plt.xscale("log")
#plt.yscale("log")

plt.show()

#plt.plot(x,lst1D[int(i*50):int(i*350)], color= "blue", label="current test")


#i = 1
#x = np.linspace(1,7,i*300)
#plt.plot(x,matrices[-1][0][int(i*50):int(i*350)], color = "red", label="RK4")
"""
i = 1
x = np.linspace(1,7,i*300)
plt.plot(x,matrices2[-1][0][int(i*50):int(i*350)], color = "green", label="RK4")
"""


#data = np.loadtxt('jan.txt')
"""
data = np.loadtxt('iV.txt')
    
Er = np.loadtxt('iEr.txt')

Ez = np.loadtxt('iEz1.txt')

r = []
z = []

dr = 1
dz = 1

epsi = 8.85e-12;

for i in range (data.shape[0]):
    r.append((0.5 + i)*dr)

#for j in range(data.shape[1]):
#    z.append((0.5 + j)*dz)

def jan(r,z):
    return math.exp(-r*r -z*z)

def ErJAN(r,z):
    return 2 * r * jan(r,z)

def EzJAN(r,z):
    return 2 * z * jan(r,z)

teste = []

Er_teste = []
Ez_teste = []

for i in r:

    row = []
    rowEr = []
    rowEz = []
    for j in z:

        row.append(jan(i,j))
        rowEr.append(ErJAN(i,j))
        rowEz.append(EzJAN(i,j))

#    rowEz.pop()
    teste.append(row)
    Er_teste.append(rowEr)
    Ez_teste.append(rowEz)

Er_teste.pop()

matrix = np.array(teste)
matrixEr = np.array(Er_teste)
matrixEz = np.array(Ez_teste)



# Calculate the global min and max values across all plots
#vmin = min(data.min(), matrix.min())
#vmax = max(data.max(), matrix.max())

#vminEr = min(Er.min(), matrixEr.min())
#vmaxEr = max(Er.max(), matrixEr.max())

#vminEz = min(Ez.min(), matrixEz.min())
#vmaxEz = max(Ez.max(), matrixEz.max())

plt.figure(figsize=(15, 5))

# Plot 1: Density Plot
plt.subplot(2, 3, 1)
scatter = plt.imshow(data, cmap='pink', interpolation='nearest')#, vmin=vmin, vmax=vmax)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
#plt.title(r'Simulation for $\rho = (6 - 4z^2 - 4r^2)e^{-z^2-r^2} $')
plt.title(r'Simulation for no charge density')
plt.gca().invert_yaxis()

#print (data -matrix)
# Show the plot
# Plot 1: Density Plot
plt.subplot(2, 3, 2)
#scatter = plt.imshow(Er, cmap='pink', interpolation='nearest')#, vmin=vminEr, vmax=vmaxEr)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
plt.title(r'Er simul')
plt.gca().invert_yaxis()

plt.subplot(2, 3, 3)
#scatter = plt.imshow(Ez, cmap='pink', interpolation='nearest')#, vmin=vminEz, vmax=vmaxEz)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
plt.title(r'Ez simul')
plt.gca().invert_yaxis()

#data = np.loadtxt('jan.txt')
data = np.loadtxt('fV.txt')

Er = np.loadtxt('fEr.txt')
"""
Ez = np.loadtxt('fEz1.txt')
"""
print(Ez[0])

# Plot 1: Density Plot
plt.subplot(2, 3, 4)
scatter = plt.imshow(data, cmap='pink', interpolation='nearest')#, vmin=vmin, vmax=vmax)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
#plt.title(r'Simulation for $\rho = (6 - 4z^2 - 4r^2)e^{-z^2-r^2} $')
plt.title(r'Simulation for no charge density')
plt.gca().invert_yaxis()


# Show the plot
# Plot 1: Density Plot
plt.subplot(2, 3, 5)
#scatter = plt.imshow(Er, cmap='pink', interpolation='nearest')#, vmin=vminEr, vmax=vmaxEr)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
plt.title(r'Er simul')
plt.gca().invert_yaxis()

plt.subplot(2, 3, 6)
#scatter = plt.imshow(Ez, cmap='pink', interpolation='nearest')#, vmin=vminEz, vmax=vmaxEz)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
plt.title(r'Ez simul')
plt.gca().invert_yaxis()

plt.show()
"""
def load_all_matrices(filename):
    nmax = -1e50;
    nmin = 1e50;
    matrices = []
    times = []
    with open(filename, 'r') as f:
        matrix_data = []
        for line in f:
            if line.startswith("Time:"):
                times.append(float(line.split()[1]))
                if matrix_data:    
                    matrices.append(np.array(matrix_data))
                    matrix_data = []
            elif line.strip() == "----":
                continue
            else:
                row = [float(val) for val in line.split()]

                if (max(row) > nmax):
                        nmax = max(row)
                if (min(row) < nmin):
                    nmin = min(row)

                matrix_data.append(row)
        if matrix_data:  # Add the last matrixb
            matrices.append(np.array(matrix_data))
    return times, matrices, nmax, nmin

# Load the data
times, matrices, nmax, nmin = load_all_matrices("build/rho_data.txt")
#times2, matrices2, nmax2, nmin2 = load_all_matrices("build/rho_data RK4 O.1PS NO CURRLIM UNO3.txt")
#times3, matrices3, nmax3, nmin3 = load_all_matrices("build/rho_data UNO3 norm.txt")
#times4, matrices4, nmax4, nmin4 = load_all_matrices("build/rho_data UNO3 trap.txt")
#times5, matrices5, nmax5, nmin5 = load_all_matrices("build/rho_data UNO3 RK4.txt")
#
#times6, matrices6, nmax6, nmin6 = load_all_matrices("build/rho_data SOLUTION SP.txt")
#times7, matrices7, nmax7, nmin7 = load_all_matrices("build/rho_data SP norm.txt")
#times8, matrices8, nmax8, nmin8 = load_all_matrices("build/rho_data SP trap.txt")
#times9, matrices9, nmax9, nmin9 = load_all_matrices("build/rho_data SP RK4.txt")

# Set up the figure and axis
#fig, ax = plt.subplots()
#im = ax.imshow(matrices[0],vmin=nmin, vmax=nmax)# cmap='viridis',norm=SymLogNorm(linthresh=1, vmin=-nmin, vmax=nmax))  # vmin=nmin, vmax=nmax,Adjust vmin/vmax as needed
#
#def update(frame):
#    im.set_array(matrices[frame])
#    ax.set_title(f'Time = {times[frame]}')
#    return [im]
#
## Create an animation
#ani = animation.FuncAnimation(fig, update, frames=len(matrices), blit=True)
#plt.colorbar(im)
#plt.show()
"""
# Define the two new lists
list1 = [
    -1.58181e-23, -1.85368e-24, -2.17228e-25, -2.54564e-26, -2.98317e-27,
    -3.49591e-28, -4.09677e-29, -4.8009e-30, -5.62605e-31, -6.59303e-32,
    -7.72621e-33, -9.05415e-34, -1.06103e-34, -1.2434e-35, -1.45711e-36,
    -1.70755e-37, -2.00103e-38, -2.34496e-39, -2.748e-40, -3.22031e-41,
    -3.7738e-42, -4.42243e-43, -5.18253e-44, -1.14348e-44, -1.14348e-44,
    -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6,
    -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6,
    -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6,
    -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6,
    -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6,
    -25.5999, -25.5994, -25.5952, -25.5588, -25.2484, -22.6, -3, -0.351562,
    -0.0411987, -0.00482798, -0.000565778, -6.63022e-05, -7.76978e-06,
    -9.10522e-07, -1.06702e-07, -1.25041e-08, -1.46533e-09, -1.71718e-10,
    -2.01232e-11, -2.35819e-12, -2.7635e-13, -3.23848e-14, -3.79509e-15,
    -4.44737e-16, -5.21176e-17, -6.10753e-18, -7.15726e-19, -8.38742e-20,
    -9.82901e-21, -1.15184e-21, -1.34981e-22
]

list2 = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6,
    -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6,
    -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6,
    -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6, -25.6,
    -25.6, -25.6, -25.6, -25.6, -25.6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
]

x = []
for i in range(100):
    x.append(i)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(np.array(list1) * -1, label='After', color='blue')
plt.plot(np.array(list2) * -1, label='Before', color='orange')

# Adding titles and labels
plt.title("Overlay of Two Lists")
plt.xlabel("Index")
plt.ylabel("Value")
plt.legend()

# Show the plot
plt.show()
"""
"""
# Define the two new lists
list1 = [
    -8.38742e-021, -9.82901e-022, -1.15184e-022, -1.34981e-023, s
    -1.58181e-024, -1.85368e-025, -2.17228e-026, -2.54564e-027, 
    -2.98317e-028, -3.49591e-029, -4.09677e-030, -4.8009e-031, 
    -5.62605e-032, -6.59303e-033, -7.72621e-034, -9.05415e-035, 
    -1.06103e-035, -1.2434e-036, -1.45711e-037, -1.70755e-038, 
    -1.70755e-038, -2.56, -5.12, -7.68, -10.24, -12.8, -15.36, 
    -17.92, -20.48, -23.04, -25.6, -28.16, -30.72, -33.28, -35.84, 
    -38.4, -40.96, -43.52, -46.08, -48.64, -51.2, -53.76, -56.32, 
    -58.88, -61.44, -64, -66.56, -69.12, -71.68, -73.94, -73.94, 
    -71.68, -69.12, -66.56, -64, -61.44, -58.88, -56.32, -53.76, 
    -51.2, -48.64, -46.08, -43.52, -40.96, -38.4, -35.84, -33.28, 
    -30.72, -28.16, -25.6, -23.04, -20.48, -17.92, -15.36, -12.8, 
    -10.24, -7.68, -5.12, -2.56, -0.3, -0.0351562, -0.00411987, 
    -0.000482798, -5.65778e-005, -6.63022e-006, -7.76978e-007, 
    -9.10522e-008, -1.06702e-008, -1.25041e-009, -1.46533e-010, 
    -1.71718e-011, -2.01232e-012, -2.35819e-013, -2.7635e-014, 
    -3.23848e-015, -3.79509e-016, -4.44737e-017, -5.21176e-018, 
    -6.10753e-019, -7.15726e-020
]


list2 = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    -2.56, -5.12, -7.68, -10.24, -12.8, -15.36, -17.92, -20.48, -23.04, 
    -25.6, -28.16, -30.72, -33.28, -35.84, -38.4, -40.96, -43.52, -46.08, 
    -48.64, -51.2, -53.76, -56.32, -58.88, -61.44, -64, -66.56, -69.12, 
    -71.68, -74.24, -74.24, -71.68, -69.12, -66.56, -64, -61.44, -58.88, 
    -56.32, -53.76, -51.2, -48.64, -46.08, -43.52, -40.96, -38.4, -35.84, 
    -33.28, -30.72, -28.16, -25.6, -23.04, -20.48, -17.92, -15.36, -12.8, 
    -10.24, -7.68, -5.12, -2.56, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0
]


x = []
for i in range(100):
    x.append(i)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(np.array(list1) * -1, label='After', color='blue')
plt.plot(np.array(list2) * -1, label='Before', color='orange')

# Adding titles and labels
plt.title("Overlay of Two Lists")
plt.xlabel("Index")
plt.ylabel("Value")
plt.legend()

# Show the plot
plt.show()
"""
"""
list1 = [
    -160, -160, -160, -160, -160, -160, -160, -160, -160, -160, 
-160, -160, -160, -160, -160, -160, -160, -160, -160, -160, 
-160, -160, -160, -160, -160, -160, -160, -160, -160, -160, 
-160, -160, -160, -160, -160, -160, -160, -160, -160, -160, 
-160, -160, -160, -160, -160, -160, -160, -160, -160, -160, -160]

list2 = [-3053.53, -39.6823, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257,
        -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257,
        -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257,
        -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257,
        -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -39.4257, -3053.53]

list3 = [-3976.09, -14.7546, -0.51802, -0.39366, -0.392561, -0.392551, -0.392551, -0.392551, -0.392551, -0.392551,
        -0.392551, -0.392551, -0.392551, -0.392551, -0.392551, -0.392551, -0.392551, -0.392551, -0.392551, -0.392551,
        -0.392551, -0.392551, -0.392551, -0.392551, -0.392551, -0.392551, -0.392551, -0.392551, -0.392551, -0.392551,
        -0.392551, -0.392551, -0.392551, -0.392551, -0.392551, -0.392551, -0.392551, -0.392551, -0.392551, -0.392551,
        -0.392551, -0.392551, -0.392551, -0.392551, -0.392561, -0.392561,-0.39366, -0.51802, -14.7546, -3976.09]


sum1 = 0
sum2 = 0
"""
"""
for i in range(len(list1)):
    sum1 += list1[i]
    sum2 += list2[i]
"""

"""
print(sum1)
print(sum2)

x = []
for i in range(50):
    x.append(i)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(np.array(list1) * -1, label='After', color='blue')
plt.plot(np.array(list2) * -1, label='Middle', color='orange')
plt.plot(np.array(list3) * -1, label='at 2.24829e-008s', color='red')

# Adding titles and labels
plt.title("Overlay of Two Lists")
plt.xlabel("Index")
plt.ylabel("Value")
plt.legend()
 
plt.yscale("log")
plt.ylim(0,5e3)
# Show the plot
plt.show()
"""
"""
dt_data = np.loadtxt("build/dt_data.txt")

plt.plot(dt_data, linestyle='-', label='Data', lw=2)
plt.title("Timestep Per iteration")
plt.xlabel("Iteration")
plt.ylabel("dt(s)")
plt.legend()
plt.grid(True)
plt.show()
"""
i = 1
x = np.linspace(1,7,int(i*300))
x_original = np.linspace(0,10,int(500))
# FOR UNO3
plt.plot(x,np.array(matrices[-1][-1][int(i*50):int(i*350)]), color= "black", label="current test", ls = "solid")
#plt.plot(x,np.array(matrices[-1][-2][int(i*50):int(i*350)]), color= "red", label="current test", ls = "--")

#plt.plot(x_original,np.array(matrices2[-1][0][int(0):int(500)]), color= "green", label="sol UNO3", lw= 1)
#plt.plot(x_original,np.array(matrices3[-1][0][int(i*50):int(i*350)]), color= "blue", label="nor", lw= 1)
#plt.plot(x_original,np.array(matrices4[-1][0][int(i*50):int(i*350)]), color= "red", label="trap", lw= 1)
#plt.plot(x_original,np.array(matrices5[-1][0][int(i*50):int(i*350)]), color= "black", label="RK4", lw= 1, ls= "--")

# for SUPERBEE
#plt.plot(x_original,np.array(matrices6[-1][0][int(i*50):int(i*350)]), color= "black", label="sol SUPERBEE", lw= 1)
#plt.plot(x_original,np.array(matrices7[-1][0][int(i*50):int(i*350)]), color= "blue", label="nor", lw= 1)
#plt.plot(x_original,np.array(matrices8[-1][0][int(i*50):int(i*350)]), color= "red", label="trap", lw= 1)
#plt.plot(x_original,np.array(matrices9[-1][0][int(i*50):int(i*350)]), color= "black", label="RK4", lw= 1, ls= "--")
plt.legend()
x = np.linspace(0,10,int(i*500))
#plt.plot(x,np.array(matrices[-2][0][int(i*50):int(i*350)]))
#plt.plot(x,np.array(matrices[-3][0][int(i*50):int(i*350)]))
#plt.plot(x,np.array(matrices[-4][0][int(i*50):int(i*350)]))
#plt.plot(x,np.array(matrices[-5][0][int(i*50):int(i*350)]))
#plt.plot(x,np.array(matrices[-6][0][int(i*50):int(i*350)]))
print(matrices[-1][0])

i = 1
x = np.linspace(1,7,int(i*500))
#plt.plot(x,matrices[-1][-1][int(i*50):int(i*350)], color = "red", label="RK4")
"""
i = 1
x = np.linspace(1,7,i*300)
plt.plot(x,matrices2[-1][0][int(i*50):int(i*350)], color = "green", label="RK4")
"""
plt.yscale("log")
plt.ylim(1e12,1e21)
plt.grid(True)
plt.title(r'Density(z)')
plt.show()

plt.title(r'E(z)')
plt.grid(True)

plt.plot(x,Ez[int(i*0):int(i*500)])
#plt.yscale("log")
#splt.ylim(1e1,1e21)
plt.show()