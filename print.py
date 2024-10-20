import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
from matplotlib.colors import SymLogNorm

#data = np.loadtxt('jan.txt')
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

for j in range(data.shape[1]):
    z.append((0.5 + j)*dz)

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

    rowEz.pop()
    teste.append(row)
    Er_teste.append(rowEr)
    Ez_teste.append(rowEz)

Er_teste.pop()

matrix = np.array(teste)
matrixEr = np.array(Er_teste)
matrixEz = np.array(Ez_teste)

# Calculate the global min and max values across all plots
vmin = min(data.min(), matrix.min())
vmax = max(data.max(), matrix.max())

vminEr = min(Er.min(), matrixEr.min())
vmaxEr = max(Er.max(), matrixEr.max())

vminEz = min(Ez.min(), matrixEz.min())
vmaxEz = max(Ez.max(), matrixEz.max())

plt.figure(figsize=(15, 5))

# Plot 1: Density Plot
plt.subplot(1, 3, 1)
scatter = plt.imshow(data, cmap='pink', interpolation='nearest')#, vmin=vmin, vmax=vmax)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
#plt.title(r'Simulation for $\rho = (6 - 4z^2 - 4r^2)e^{-z^2-r^2} $')
plt.title(r'Simulation for no charge density')
plt.gca().invert_yaxis()

print (data -matrix)
# Show the plot
# Plot 1: Density Plot
plt.subplot(1, 3, 2)
scatter = plt.imshow(Er, cmap='pink', interpolation='nearest')#, vmin=vminEr, vmax=vmaxEr)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
plt.title(r'Er simul')
plt.gca().invert_yaxis()

plt.subplot(1, 3, 3)
scatter = plt.imshow(Ez, cmap='pink', interpolation='nearest')#, vmin=vminEz, vmax=vmaxEz)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
plt.title(r'Ez simul')
plt.gca().invert_yaxis()

plt.show()

def load_all_matrices(filename):
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
                matrix_data.append(row)
        if matrix_data:  # Add the last matrixb
            matrices.append(np.array(matrix_data))
    return times, matrices

# Load the data
times, matrices = load_all_matrices("build/rho_data.txt")

# Set up the figure and axis
fig, ax = plt.subplots()
im = ax.imshow(matrices[0], cmap='viridis')#,norm=SymLogNorm(linthresh=1))  # Adjust vmin/vmax as needed

def update(frame):
    im.set_array(matrices[frame])
    ax.set_title(f'Time = {times[frame]}')
    return [im]

# Create an animation
ani = animation.FuncAnimation(fig, update, frames=len(matrices), blit=True)
plt.colorbar(im)
plt.show()

#data = np.loadtxt('jan.txt')
data = np.loadtxt('fV.txt')

Er = np.loadtxt('fEr.txt')

Ez = np.loadtxt('fEz1.txt')

r = []
z = []

dr = 1
dz = 1

epsi = 8.85e-12;

for i in range (data.shape[0]):
    r.append((0.5 + i)*dr)

for j in range(data.shape[1]):
    z.append((0.5 + j)*dz)

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

    rowEz.pop()
    teste.append(row)
    Er_teste.append(rowEr)
    Ez_teste.append(rowEz)

Er_teste.pop()

matrix = np.array(teste)
matrixEr = np.array(Er_teste)
matrixEz = np.array(Ez_teste)

# Calculate the global min and max values across all plots
vmin = min(data.min(), matrix.min())
vmax = max(data.max(), matrix.max())

vminEr = min(Er.min(), matrixEr.min())
vmaxEr = max(Er.max(), matrixEr.max())

vminEz = min(Ez.min(), matrixEz.min())
vmaxEz = max(Ez.max(), matrixEz.max())

plt.figure(figsize=(15, 5))

# Plot 1: Density Plot
plt.subplot(1, 3, 1)
scatter = plt.imshow(data, cmap='pink', interpolation='nearest')#, vmin=vmin, vmax=vmax)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
#plt.title(r'Simulation for $\rho = (6 - 4z^2 - 4r^2)e^{-z^2-r^2} $')
plt.title(r'Simulation for no charge density')
plt.gca().invert_yaxis()

print (data -matrix)
# Show the plot
# Plot 1: Density Plot
plt.subplot(1, 3, 2)
scatter = plt.imshow(Er, cmap='pink', interpolation='nearest')#, vmin=vminEr, vmax=vmaxEr)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
plt.title(r'Er simul')
plt.gca().invert_yaxis()

plt.subplot(1, 3, 3)
scatter = plt.imshow(Ez, cmap='pink', interpolation='nearest')#, vmin=vminEz, vmax=vmaxEz)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
plt.title(r'Ez simul')
plt.gca().invert_yaxis()


plt.show()


def load_all_matrices(filename):
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
                matrix_data.append(row)
        if matrix_data:  # Add the last matrixb
            matrices.append(np.array(matrix_data))
    return times, matrices

# Load the data
times, matrices = load_all_matrices("build/rho_data.txt")

# Set up the figure and axis
fig, ax = plt.subplots()
im = ax.imshow(matrices[0], cmap='viridis')#,norm=SymLogNorm(linthresh=1))  # Adjust vmin/vmax as needed

def update(frame):
    im.set_array(matrices[frame])
    ax.set_title(f'Time = {times[frame]}')
    return [im]

# Create an animation
ani = animation.FuncAnimation(fig, update, frames=len(matrices), blit=True)
plt.colorbar(im)
plt.show()