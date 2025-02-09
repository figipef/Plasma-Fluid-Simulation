import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
from matplotlib.colors import SymLogNorm
from scipy.optimize import curve_fit

# Initialize lists to store the columns
column1 = []
column2 = []

# Open and read the txt file
with open("bolsig/N2mob.txt", "r") as file:
    for line in file:
        # Split each line into two parts
        values = line.split()
        # Convert the values to floats and store them in the respective lists
        column1.append(float(values[0]))
        column2.append(float(values[1]))

x = np.logspace(-3, 4,100)

def arctan(x,a,b,c,d):
    return a * np.arctan(np.log10(x) * b - c) + d

# Fit the data with free parameters
popt, pcov = curve_fit(
    arctan, 
    column1, 
    column2, 
    p0=[-1e25, 1, -5,1e25]  # Initial guesses for a, b, and c
)

# Extract fitted parameters
a_fit, b_fit, c_fit, d_fit = popt
print(a_fit, b_fit, c_fit, d_fit)

#print(f)
plt.plot(x, arctan(x, a_fit, b_fit, c_fit, d_fit), 
         label=f'Fitted arctan: a={a_fit:.2e}, b={b_fit:.2e}, c={c_fit:.2e}', 
         color='red')
plt.scatter(column1, column2)
plt.plot()
plt.xscale("log")
plt.yscale("log")
plt.show()

coeffs = [ 4.19038e+024,
-2.50599e+024,
 3.44633e+023]

# Initialize lists to store the columns
column1 = []
column2 = []

# Open and read the txt file
with open("bolsig/N2temp.txt", "r") as file:
    for line in file:
        # Split each line into two parts
        values = line.split()
        # Convert the values to floats and store them in the respective lists
        column1.append(float(values[0]))
        column2.append(float(values[1]))

#print(column1[0])
x = np.logspace(-4, 4,1000)

def dif(x):

    if x < 1:
        return 1e24
    else:
        return 1e24 * (np.exp(0.02*np.log( (x*100))**2) - 0.5)

def difusion(x,a1,a2,b2,c2,d2,a3,b3,c3,d3,e3,a4,b4,c4):
    return np.piecewise(x,
        [x <= 0.1, (x >0.1) & (x <= 10*np.sqrt(10)), (x > 10 *np.sqrt(10)) & (x <=100), x > 100],
        [
            lambda x: a1,
            lambda x: a2*np.arctan(b2*np.log10(x) - c2) + d2,
            lambda x: a3*(np.log10(x) -b3)**2 + c3 * (np.log10(x)-d3) + e3,
            lambda x: a4* np.exp(b4 * np.log10(x)) + c4
        ])

p0=[9.798e23, 2.9e23, 2.818,2.265e-1,1.397e24,3.770921e24,1.66159,-100,1,1.7099e24,4.918e24,3.711e-1,-8.059e24]  # Initial guesses for a, b, and c

popt, pcov = curve_fit(
    difusion, 
    column1, 
    column2, 
    p0=p0  # Initial guesses for a, b, and c
)

a1,a2,b2,c2,d2,a3,b3,c3,d3,e3,a4,b4,c4 = p0

print(a1,a2,b2,c2,d2,a3,b3,c3,d3,e3,a4,b4,c4)

plt.plot(x, difusion(x,a1,a2,b2,c2,d2,a3,b3,c3,d3,e3,a4,b4,c4), 
         label=f'Fitted arctan: a={a_fit:.2e}, b={b_fit:.2e}, c={c_fit:.2e}', 
         color='red')
plt.scatter(column1, column2)
plt.plot()
plt.xscale("log")
plt.yscale("log")
plt.xlim(1e-4, 10e3)

plt.show()

# Initialize lists to store the columns
column1 = []
column2 = []

def b(x):

    if x < 10:
        return 0
    else:
        return 1e-13 * (np.exp(-0.5*( np.log(x*0.0001))**2))

# Open and read the txt file
with open("bolsig/N2iorate.txt", "r") as file:
    for line in file:
        # Split each line into two parts
        values = line.split()
        # Convert the values to floats and store them in the respective lists
        column1.append(float(values[0]))
        column2.append(float(values[1]))

x = np.logspace(-3, 4,100)


#f = pol(x, coeffs, 3,1)
f = []
for i in x:
    #print(i)
    #print(dif(i))
    f.append(float(b(i)))

#print(f)
plt.plot(x, f)
plt.scatter(column1, column2)
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
times2, matrices2, nmax2, nmin2 = load_all_matrices("build/rho_data RK4 O.1PS NO CURRLIM UNO3.txt")
times3, matrices3, nmax3, nmin3 = load_all_matrices("build/rho_data UNO3 norm.txt")
times4, matrices4, nmax4, nmin4 = load_all_matrices("build/rho_data UNO3 trap.txt")
times5, matrices5, nmax5, nmin5 = load_all_matrices("build/rho_data UNO3 RK4.txt")

times6, matrices6, nmax6, nmin6 = load_all_matrices("build/rho_data SOLUTION SP.txt")
times7, matrices7, nmax7, nmin7 = load_all_matrices("build/rho_data SP norm.txt")
times8, matrices8, nmax8, nmin8 = load_all_matrices("build/rho_data SP trap.txt")
times9, matrices9, nmax9, nmin9 = load_all_matrices("build/rho_data SP RK4.txt")

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

dt_data = np.loadtxt("build/dt_data.txt")

plt.plot(dt_data, linestyle='-', label='Data', lw=2)
plt.title("Timestep Per iteration")
plt.xlabel("Iteration")
plt.ylabel("dt(s)")
plt.legend()
plt.grid(True)
plt.show()

i = 1
x = np.linspace(0,10,int(i*500))
x_original = np.linspace(0,10,int(500))
# FOR UNO3
plt.plot(x,np.array(matrices[-1][-1][int(i*0):int(i*500)]), color= "black", label="current test", ls = "solid")
#plt.plot(x,np.array(matrices[-1][-2][int(i*50):int(i*350)]), color= "red", label="current test", ls = "--")

plt.plot(x_original,np.array(matrices2[-1][0][int(0):int(500)]), color= "green", label="sol UNO3", lw= 1)
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
x = np.linspace(0,10,int(i*500))
#plt.plot(x,matrices[-1][-1][int(i*50):int(i*350)], color = "red", label="RK4")
"""
i = 1
x = np.linspace(1,7,i*300)
plt.plot(x,matrices2[-1][0][int(i*50):int(i*350)], color = "green", label="RK4")
"""
plt.yscale("log")
plt.ylim(1e17,1e21)
plt.grid(True)
plt.title(r'Density(z)')
plt.show()

plt.title(r'E(z)')
plt.grid(True)

plt.plot(x,Ez[int(i*0):int(i*500)])
#plt.yscale("log")
#splt.ylim(1e1,1e21)
plt.show()