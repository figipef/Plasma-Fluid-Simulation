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
plt.subplot(2, 3, 1)
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
plt.subplot(2, 3, 2)
scatter = plt.imshow(Er, cmap='pink', interpolation='nearest')#, vmin=vminEr, vmax=vmaxEr)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
plt.title(r'Er simul')
plt.gca().invert_yaxis()

plt.subplot(2, 3, 3)
scatter = plt.imshow(Ez, cmap='pink', interpolation='nearest')#, vmin=vminEz, vmax=vmaxEz)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
plt.title(r'Ez simul')
plt.gca().invert_yaxis()

#data = np.loadtxt('jan.txt')
data = np.loadtxt('fV.txt')

Er = np.loadtxt('fEr.txt')

Ez = np.loadtxt('fEz1.txt')

# Plot 1: Density Plot
plt.subplot(2, 3, 4)
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
plt.subplot(2, 3, 5)
scatter = plt.imshow(Er, cmap='pink', interpolation='nearest')#, vmin=vminEr, vmax=vmaxEr)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
plt.title(r'Er simul')
plt.gca().invert_yaxis()

plt.subplot(2, 3, 6)
scatter = plt.imshow(Ez, cmap='pink', interpolation='nearest')#, vmin=vminEz, vmax=vmaxEz)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
plt.title(r'Ez simul')
plt.gca().invert_yaxis()

plt.show()

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

# Set up the figure and axis
fig, ax = plt.subplots()
im = ax.imshow(matrices[0],vmin=nmin, vmax=nmax)# cmap='viridis',norm=SymLogNorm(linthresh=1, vmin=-nmin, vmax=nmax))  # vmin=nmin, vmax=nmax,Adjust vmin/vmax as needed

def update(frame):
    im.set_array(matrices[frame])
    ax.set_title(f'Time = {times[frame]}')
    return [im]

# Create an animation
ani = animation.FuncAnimation(fig, update, frames=len(matrices), blit=True)
plt.colorbar(im)
plt.show()
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
    -8.38742e-021, -9.82901e-022, -1.15184e-022, -1.34981e-023, 
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
list1 = [
    -160, -160, -160, -160, -160, -160, -160, -160, -160, -160, 
-160, -160, -160, -160, -160, -160, -160, -160, -160, -160, 
-160, -160, -160, -160, -160, -160, -160, -160, -160, -160, 
-160, -160, -160, -160, -160, -160, -160, -160, -160, -160, 
-160, -160, -160, -160, -160, -160, -160, -160, -160, -160, -160]

list2 = [
    -138396, -119.646, -115.763, -115.763, -115.763, -115.763, -115.763, 
    -115.763, -115.763, -115.763, -115.763, -115.763, -115.763, -115.763, 
    -115.763, -115.763, -115.763, -115.763, -115.763, -115.763, -115.763, 
    -115.763, -115.763, -115.763, -115.763, -115.763, -115.763, -115.763, 
    -115.763, -115.763, -115.763, -115.763, -115.763, -115.763, -115.763, 
    -115.763, -115.763, -115.763, -115.763, -115.763, -115.763, -115.763, 
    -115.763, -115.763, -115.763, -115.763, -115.763, -115.763, -115.763, 
    -119.646, -138396
]

sum1 = 0
sum2 = 0

for i in range(len(list1)):
    sum1 += list1[i]
    sum2 += list2[i]

print(sum1)
print(sum2)

x = []
for i in range(50):
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

plt.yscale('log') 
plt.ylim(1e0,1e6)
# Show the plot
plt.show()

list1 = [
    -160, -160, -160, -160, -160, -160, -160, -160, -160, -160, 
-160, -160, -160, -160, -160, -160, -160, -160, -160, -160]

list2 = [
    -1598.42, 1154.01, -1153.93, 802.903, -803.859, 389.576, -381.919, 73.2396,
    -61.2774, -20.7472, -20.7472, -61.2774, 73.2396, -381.919, 389.576, 
    -803.856, 802.894, -1154, 1154.26, -1598.43
]

sum1 = 0
sum2 = 0

for i in range(len(list1)):
    sum1 += list1[i]
    sum2 += list2[i]

print(sum1)
print(sum2)

x = []
for i in range(50):
    x.append(i)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(np.array(list1) * -1, label='Before', color='blue')
plt.plot(np.array(list2) * -1, label='After', color='orange')

# Adding titles and labels
plt.title("Overlay of Two Lists")
plt.xlabel("Index")
plt.ylabel("Value")
plt.legend()

# Show the plot
plt.show()