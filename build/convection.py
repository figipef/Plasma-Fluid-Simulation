import matplotlib.pyplot as plt
import math

def create_square_function_list(size, a):
    # Initialize list with zeros
    square_list = [0] * size
    
    # Define the range for non-zero values
    start = size // 4
    end = 3 * size // 4
    
    # Set values to `a` within the specified range
    for i in range(150, 350):
        square_list[i] = a
    """
    for i in range(300, 350):
        square_list[i] = a*(i-300)/(50)
        square_list[199-i] = a*(i-300)/(50)
"""
    return square_list

def sign(a,b):
    if a>b:
        return 1
    elif b>a:
        return -1
    else:
        return 0

# Input for list size and the value 'a'
size = 500
a = 1

# Create the square function list
square_list = create_square_function_list(size, a)

plt.plot(square_list, label="Square Function")

new_square_list = square_list


n = 3195

dx = 1
mu = 0.625

u = 1

dt = mu * dx/u

last_phi_E = 0
"""
for i in range(n):
    first_value = new_square_list[0]
    last_value = new_square_list[size-1]

    
    for j in range(size):

        if j == 0:
            phi_W = last_phi_E
            phi_E = new_square_list[j]  + 0.5 * sign(new_square_list[j+1], new_square_list[j]) * (dx - u * dt)*min(math.fabs((new_square_list[j] - last_value)/dx), math.fabs((new_square_list[j+1] - new_square_list[j])/dx))

        if j == size - 1:
            phi_W = last_phi_E
            phi_E = new_square_list[j]  + 0.5 * sign(first_value, new_square_list[j]) * (dx - u * dt)*min(math.fabs((new_square_list[j] - new_square_list[j-1])/dx), math.fabs((first_value - new_square_list[j])/dx))

        if j > 0 and j < size - 1:

            #if min(math.fabs((new_square_list[j] - new_square_list[j-1])/dx), math.fabs((new_square_list[j+1] - new_square_list[j])/dx)) != 0:
            #print(math.fabs((new_square_list[j] - new_square_list[j-1])/dx),math.fabs((new_square_list[j+1] - new_square_list[j])/dx) )    
            #print(min(math.fabs((new_square_list[j] - new_square_list[j-1])/dx), math.fabs((new_square_list[j+1] - new_square_list[j])/dx)))
            phi_W = last_phi_E
            phi_E = new_square_list[j]  + 0.5 * sign(new_square_list[j+1], new_square_list[j]) * (dx - u * dt)*min(math.fabs((new_square_list[j] - new_square_list[j-1])/dx), math.fabs((new_square_list[j+1] - new_square_list[j])/dx))

        #print(phi_W, phi_E)
        new_square_list[j] = new_square_list[j] + (u*phi_W - u*phi_E) * dt/dx
        last_phi_E = phi_E
"""

for i in range(n):
    # Keep a copy of the current state to avoid modifying while iterating
    temp_list = new_square_list.copy()
    
    # Loop through the grid
    for j in range(size):
        if j == 0:
            #phi_W = temp_list[-1]
            phi_W = temp_list[-1] + sign(temp_list[j], temp_list[- 1]) * (dx - u * dt) * abs((temp_list[j] - temp_list[- 1]) * (temp_list[j + 1] - temp_list[j]) / (dx*dx))/( abs((temp_list[j] - temp_list[-1]) / dx) + abs((temp_list[j + 1] - temp_list[j]) / dx) + 1e-308) # Periodic boundary on the left side
            phi_E = temp_list[j] + sign(temp_list[j + 1], temp_list[j]) * (dx - u * dt) * abs((temp_list[j] - temp_list[-1]) * (temp_list[j + 1] - temp_list[j]) / (dx*dx))/( abs((temp_list[j] - temp_list[-1]) / dx) + abs((temp_list[j + 1] - temp_list[j]) / dx) + 1e-308) #min(abs((temp_list[j] - temp_list[-1]) / dx), abs((temp_list[j + 1] - temp_list[j]) / dx))
            #phi_E = temp_list[j] + 0.5 * sign(temp_list[j + 1], temp_list[j]) * (dx - u * dt) * min(abs((temp_list[j] - temp_list[-1]) / dx), abs((temp_list[j + 1] - temp_list[j]) / dx))


        elif j == size - 1:
            phi_W = temp_list[j - 1] + sign(temp_list[j], temp_list[j - 1]) * (dx - u * dt)* abs((temp_list[j] - temp_list[j - 1]) * (temp_list[0] - temp_list[j]) / (dx*dx))/( abs((temp_list[j] - temp_list[j-1]) / dx) + abs((temp_list[0] - temp_list[j]) / dx) + 1e-308) # * min(abs((temp_list[j] - temp_list[j - 1]) / dx), abs((temp_list[0] - temp_list[j]) / dx))
            #phi_W = temp_list[j - 1] + 0.5 * sign(temp_list[j], temp_list[j - 1]) * min(abs((temp_list[j] - temp_list[j - 1]) / dx), abs((temp_list[0] - temp_list[j]) / dx))
            #phi_E = temp_list[0]  # Periodic boundary on the right side
            phi_E = temp_list[j] + sign(temp_list[0], temp_list[j]) * (dx - u * dt) * abs((temp_list[j] - temp_list[j - 1]) * (temp_list[0] - temp_list[j]) / (dx*dx))/( abs((temp_list[j] - temp_list[j-1]) / dx) + abs((temp_list[0] - temp_list[j]) / dx) + 1e-308)#min(abs((temp_list[j] - temp_list[j - 1]) / dx), abs((temp_list[j + 1] - temp_list[j]) / dx))

        else:
            phi_W = temp_list[j - 1] + sign(temp_list[j], temp_list[j - 1]) * (dx - u * dt) * abs((temp_list[j] - temp_list[j - 1]) * (temp_list[j + 1] - temp_list[j]) / (dx*dx))/( abs((temp_list[j] - temp_list[j-1]) / dx) + abs((temp_list[j + 1] - temp_list[j]) / dx) + 1e-308) #* min(abs((temp_list[j] - temp_list[j - 1]) / dx), abs((temp_list[j + 1] - temp_list[j]) / dx))
            phi_E = temp_list[j] + sign(temp_list[j + 1], temp_list[j]) * (dx - u * dt) * abs((temp_list[j] - temp_list[j - 1]) * (temp_list[j + 1] - temp_list[j]) / (dx*dx))/( abs((temp_list[j] - temp_list[j-1]) / dx) + abs((temp_list[j + 1] - temp_list[j]) / dx) + 1e-308)#min(abs((temp_list[j] - temp_list[j - 1]) / dx), abs((temp_list[j + 1] - temp_list[j]) / dx))

            #phi_W = temp_list[j - 1] + 0.5 * sign(temp_list[j], temp_list[j - 1]) * (dx - u * dt) * min(abs((temp_list[j] - temp_list[j - 1]) / dx), abs((temp_list[j + 1] - temp_list[j]) / dx))
            #phi_E = temp_list[j] + 0.5 * sign(temp_list[j + 1], temp_list[j]) * (dx - u * dt) * min(abs((temp_list[j] - temp_list[j - 1]) / dx), abs((temp_list[j + 1] - temp_list[j]) / dx))

        # Update new_square_list using phi_W and phi_E
        new_square_list[j] = temp_list[j] + (u * phi_W - u * phi_E) * dt / dx

print(new_square_list)
# Plot the list

plt.plot(new_square_list, label="new_Square Function")
plt.title("Square Function Plot")
plt.xlabel("Index")
plt.ylabel("Value")
plt.legend()
plt.show()
"""
import numpy as np
import matplotlib.pyplot as plt

# Parameters
nx = size             # Number of spatial points
L = 1.0              # Length of domain
dx = L / (nx - 1)    # Spatial step size
v = 1.0              # Velocity of advection
CFL = 0.625           # CFL number
dt = CFL * dx / v    # Time step based on CFL

# Grid and initial condition
x = np.linspace(0, L, nx)
#u = np.zeros(nx)
u = square_list
# Initial condition: square function
#u[int(nx / 4):int(3 * nx / 4)] = 1

# UNO2 function for flux-limited advection
def uno2_flux_limiter(u, dx):
    flux = np.zeros_like(u)
    # Compute fluxes for each point except boundaries
    for i in range(1, nx - 1):
        # Left-biased and right-biased differences
        ul = u[i] - u[i - 1]
        ur = u[i + 1] - u[i]
        #print(ul, ur)
        # Calculate smoothness to determine if limiter is needed
        if ul * ur > 0:
            # Smooth region, use second-order interpolation
            flux[i] = u[i] + 0.5 * minmod(ul, ur)
        else:
            # Discontinuity detected, use first-order upwind flux
            flux[i] = u[i]


    ul = u[0] - u[nx-1]
    ur = u[1] - u[0]
   
    if ul * ur > 0:
        # Smooth region, use second-order interpolation
        flux[0] = u[0] + 0.5 * minmod(ul, ur)
    else:
        # Discontinuity detected, use first-order upwind flux
        flux[0] = u[0]

    ul = u[nx-1] - u[nx-2]
    ur = u[0] - u[nx-1]
   
    if ul * ur > 0:
        # Smooth region, use second-order interpolation
        flux[nx-1] = u[nx-1] + 0.5 * minmod(ul, ur)
    else:
        # Discontinuity detected, use first-order upwind flux
        flux[nx-1] = u[nx-1]

    return flux

# Minmod limiter to control oscillations
def minmod(a, b):
    if a * b > 0:
        return min(abs(a), abs(b)) * np.sign(a)
    else:
        return 0.0

# Time-stepping loop
nt = 3200  # Number of time steps
for _ in range(nt):
    # Compute fluxes with UNO2 scheme
    flux = uno2_flux_limiter(u, dx)
    
    # Update solution
    for i in range(1, nx):
        u[i] -= v * (flux[i] - flux[i - 1]) * dt / dx
    u[0] -= v * (flux[0] - flux[nx - 1]) * dt / dx
# Plot the result
plt.plot(u, label="UNO2 Solution")
plt.xlabel("x")
plt.ylabel("u")
plt.title("Advection using UNO2 Scheme")
plt.legend()
plt.show()
"""