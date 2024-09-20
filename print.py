import numpy as np
import matplotlib.pyplot as plt
import math

#data = np.loadtxt('jan.txt')
data = np.loadtxt('zero.txt')

Er = np.loadtxt('Er.txt')

Ez = np.loadtxt('E1z.txt')

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

# Plot 2: Correct Plot
#plt.subplot(1, 3, 2)
#scatter = plt.imshow(matrix, cmap='pink', interpolation='nearest', vmin=vmin, vmax=vmax)
#plt.colorbar(label='Data Value')
#plt.xlabel('Z_i')
#plt.ylabel('R_i')
#plt.title(r'Correct Plot $e^{-z^2-r^2}$')
#plt.gca().invert_yaxis()

# Plot 3: Difference Plot
#plt.subplot(1, 3, 3)
#scatter = plt.imshow(np.abs((data - matrix)), cmap='pink', interpolation='nearest')
#plt.colorbar(label='Data Value')
#plt.xlabel('Z_i')
#plt.ylabel('R_i')
#plt.title('Difference simul - correct')
#plt.gca().invert_yaxis()

# Show the plots
#plt.show()
#plt.figure(figsize=(15, 10))


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

# Plot 2: Correct Plot
#plt.subplot(2, 3, 2)
#scatter = plt.imshow(Er_teste, cmap='pink', interpolation='nearest', vmin=vminEr, vmax=vmaxEr)
#plt.colorbar(label='Data Value')
#plt.xlabel('Z_i')
#plt.ylabel('R_i')
#plt.title(r'Er teorico')
#plt.gca().invert_yaxis()

#plt.subplot(2, 3, 3)
#scatter = plt.imshow(np.abs((Er - Er_teste)), cmap='pink', interpolation='nearest')
#plt.colorbar(label='Data Value')
#plt.xlabel('Z_i')
#plt.ylabel('R_i')
#plt.title(r'Er difference')
#plt.gca().invert_yaxis()

plt.subplot(1, 3, 3)
scatter = plt.imshow(Ez, cmap='pink', interpolation='nearest')#, vmin=vminEz, vmax=vmaxEz)
plt.colorbar(label='Data Value')
plt.xlabel('Z_i')
plt.ylabel('R_i')
plt.title(r'Ez simul')
plt.gca().invert_yaxis()

# Plot 2: Correct Plot
#plt.subplot(2, 3, 5)
#scatter = plt.imshow(Ez_teste, cmap='pink', interpolation='nearest', vmin=vminEz, vmax=vmaxEz)
#plt.colorbar(label='Data Value')
#plt.xlabel('Z_i')
#plt.ylabel('R_i')
#plt.title(r'Ez teorico')
#plt.gca().invert_yaxis()
#
#plt.subplot(2, 3, 6)
#scatter = plt.imshow(np.abs((Ez - Ez_teste)), cmap='pink', interpolation='nearest')
#plt.colorbar(label='Data Value')
#plt.xlabel('Z_i')
#plt.ylabel('R_i')
#plt.title(r'Ez difference')
#plt.gca().invert_yaxis()

plt.show()
