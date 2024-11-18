import matplotlib.pyplot as plt
import math
import numpy as np

def create_square_function_list(size, a):
    # Initialize list with zeros
    square_list = [0] * size
    
    # Define the range for non-zero values
    start = size // 4
    end = 3 * size // 4
    
    # Set values to `a` within the specified range
    for i in range(100, 200):
        square_list[i] = a
    
    for i in range(300, 351):
        square_list[i] = a*(i-300)/(50)
        square_list[400-(i-300)] = a*(i-300)/(50)

    return square_list

def sign(a,b):
    if a>b:
        return 1
    elif b>a:
        return -1
    else:
        return 0

def mid_way_f(donor,v,dx,dt,G_C):
    return donor + 0.5 * sign(v,0)*(dx - abs(v)*dt)*G_C

# Input for list size and the value 'a'
size = 500
a = 1

# Create the square function list
square_list = create_square_function_list(size, a)

plt.plot(square_list, label="Square Function")

UnoPlusList = square_list
UnoList = square_list
Uno3MinusList = square_list
Uno3List = square_list

n = 3200

dx = 1
mu = 0.625

u = 1

dt = mu * dx/u

last_phi_E = 0

for i in range(n):
    print
    # Keep a copy of the current state to avoid modifying while iterating
    temp_list = UnoPlusList.copy()
    temp2_list = UnoList.copy()
    temp3_list = Uno3MinusList.copy()
    temp4_list = Uno3List.copy()
    # Loop through the grid
    for j in range(size):
        if j == 0:

            #UNO2 plus

            phi_W = UnoPlusList[-1] + sign(UnoPlusList[j], UnoPlusList[- 1]) * (dx - u * dt) * abs((UnoPlusList[-1] - UnoPlusList[- 2]) * (UnoPlusList[j] - UnoPlusList[-1]) / (dx*dx))/( abs((UnoPlusList[-1] - UnoPlusList[-2]) / dx) + abs((UnoPlusList[j] - UnoPlusList[-1]) / dx) + 1e-308) # Periodic boundary on the left side
            phi_E = UnoPlusList[j] + sign(UnoPlusList[j + 1], UnoPlusList[j]) * (dx - u * dt) * abs((UnoPlusList[j] - UnoPlusList[-1]) * (UnoPlusList[j + 1] - UnoPlusList[j]) / (dx*dx))/( abs((UnoPlusList[j] - UnoPlusList[-1]) / dx) + abs((UnoPlusList[j + 1] - UnoPlusList[j]) / dx) + 1e-308) #min(abs((temp_list[j] - temp_list[-1]) / dx), abs((temp_list[j + 1] - temp_list[j]) / dx))

            #UNO2

            phi_W2 = UnoList[-1] + 0.5 * sign(UnoList[j], UnoList[- 1]) * (dx - u * dt) * min(abs((UnoList[-1] - UnoList[-2]) / dx), abs((UnoList[j] - UnoList[-1]) / dx))
            phi_E2 = UnoList[j] + 0.5 * sign(UnoList[j + 1], UnoList[j]) * (dx - u * dt) * min(abs((UnoList[j] - UnoList[-1]) / dx), abs((UnoList[j + 1] - UnoList[j]) / dx))

            #UNO3 minus

            G_DC = (Uno3MinusList[j] - Uno3MinusList[-1])/dx
            G_CU = (Uno3MinusList[-1] - Uno3MinusList[-2])/dx

            if abs(G_DC - G_CU) <= 1.2 * abs(G_DC + G_CU)*0.5:
                G_C = G_DC - 2* (dx + u * dt)/3*(G_DC - G_CU)/(2*dx)
            else:
                G_C = sign(G_DC,0)*2*abs(G_DC * G_CU)/(abs(G_DC) + abs(G_CU) + 1e-308)

            phi_W3 = mid_way_f(Uno3MinusList[-1], u, dx, dt, G_C)

            G_DC = (Uno3MinusList[j+1] - Uno3MinusList[j])/dx
            G_CU = (Uno3MinusList[j] - Uno3MinusList[-1])/dx

            if abs(G_DC - G_CU) <= 1.2 * abs(G_DC + G_CU)*0.5:
                G_C = G_DC - 2* (dx + u * dt)/3*(G_DC - G_CU)/(2*dx)
            else:
                G_C = sign(G_DC,0)*2*abs(G_DC * G_CU)/(abs(G_DC) + abs(G_CU) + 1e-308)

            phi_E3 = mid_way_f(Uno3MinusList[j], u, dx, dt, G_C)

            #UNO3

            G_DC = (Uno3List[j] - Uno3List[-1])/dx
            G_CU = (Uno3List[-1] - Uno3List[-2])/dx

            if abs(G_DC - G_CU) <= 1.2 * abs(G_DC + G_CU)*0.5:
                G_C = G_DC - 2 * (dx + u * dt)/3*(G_DC - G_CU)/(2*dx)
            elif G_DC * G_CU > 0:
                G_C = sign(G_DC,0)*2*min(abs(G_DC), abs(G_CU))
            else:
                G_C = sign(G_DC,0)*min(abs(G_DC), abs(G_CU))

            phi_W4 = mid_way_f(Uno3List[-1], u, dx, dt, G_C)

            G_DC = (Uno3List[j+1] - Uno3List[j])/dx
            G_CU = (Uno3List[j] - Uno3List[-1])/dx

            if abs(G_DC - G_CU) <= 1.2 * abs(G_DC + G_CU)*0.5:
                G_C = G_DC - 2 * (dx + u * dt)/3*(G_DC - G_CU)/(2*dx)
            elif G_DC * G_CU > 0:
                G_C = sign(G_DC,0)*2*min(abs(G_DC), abs(G_CU))
            else:
                G_C = sign(G_DC,0)*min(abs(G_DC), abs(G_CU))

            phi_E4 = mid_way_f(Uno3List[j], u, dx, dt, G_C)

        elif j == size - 1:

            #UNO2 plus

            phi_W = UnoPlusList[j - 1] + sign(UnoPlusList[j], UnoPlusList[j - 1]) * (dx - u * dt)* abs((UnoPlusList[j - 1] - UnoPlusList[j - 2]) * (UnoPlusList[j] - UnoPlusList[j-1]) / (dx*dx))/( abs((UnoPlusList[j-1] - UnoPlusList[j-2]) / dx) + abs((UnoPlusList[j] - UnoPlusList[j - 1]) / dx) + 1e-308) # * min(abs((temp_list[j] - temp_list[j - 1]) / dx), abs((temp_list[0] - temp_list[j]) / dx))
            phi_E = UnoPlusList[j] + sign(UnoPlusList[0], UnoPlusList[j]) * (dx - u * dt) * abs((UnoPlusList[j] - UnoPlusList[j - 1]) * (UnoPlusList[0] - UnoPlusList[j]) / (dx*dx))/( abs((UnoPlusList[j] - UnoPlusList[j-1]) / dx) + abs((UnoPlusList[0] - UnoPlusList[j]) / dx) + 1e-308)#min(abs((temp_list[j] - temp_list[j - 1]) / dx), abs((temp_list[j + 1] - temp_list[j]) / dx))
            
            #UNO2

            phi_W2 = UnoList[j - 1] + 0.5 * sign(UnoList[j], UnoList[j - 1]) * (dx - u * dt)* min(abs((UnoList[j-1] - UnoList[j - 2]) / dx), abs((UnoList[j] - UnoList[j-1]) / dx))
            phi_E2 = UnoList[j] + 0.5 * sign(UnoList[0], UnoList[j]) * (dx - u * dt) * min(abs((UnoList[j] - UnoList[j - 1]) / dx), abs((UnoList[0] - UnoList[j]) / dx))

            #UNO3 minus

            G_DC = (Uno3MinusList[j] - Uno3MinusList[j-1])/dx
            G_CU = (Uno3MinusList[j-1] - Uno3MinusList[j-2])/dx

            if abs(G_DC - G_CU) <= 1.2 * abs(G_DC + G_CU)*0.5:
                G_C = G_DC - 2* (dx + u * dt)/3*(G_DC - G_CU)/(2*dx)
            else:
                G_C = sign(G_DC,0)*2*abs(G_DC * G_CU)/(abs(G_DC) + abs(G_CU) + 1e-308)

            phi_W3 = mid_way_f(Uno3MinusList[j-1], u, dx, dt, G_C)

            G_DC = (Uno3MinusList[0] - Uno3MinusList[j])/dx
            G_CU = (Uno3MinusList[j] - Uno3MinusList[j-1])/dx

            if abs(G_DC - G_CU) <= 1.2 * abs(G_DC + G_CU)*0.5:
                G_C = G_DC - 2* (dx + u * dt)/3*(G_DC - G_CU)/(2*dx)
            else:
                G_C = sign(G_DC,0)*2*abs(G_DC * G_CU)/(abs(G_DC) + abs(G_CU) + 1e-308)

            phi_E3 = mid_way_f(Uno3MinusList[j], u, dx, dt, G_C)

            #UNO3
            G_DC = (Uno3List[j] - Uno3List[j-1])/dx
            G_CU = (Uno3List[j-1] - Uno3List[j-2])/dx

            if abs(G_DC - G_CU) <= 1.2 * abs(G_DC + G_CU)*0.5:
                G_C = G_DC - 2 * (dx + u * dt)/3*(G_DC - G_CU)/(2*dx)
            elif G_DC * G_CU > 0:
                G_C = sign(G_DC,0)*2*min(abs(G_DC), abs(G_CU))
            else:
                G_C = sign(G_DC,0)*min(abs(G_DC), abs(G_CU))

            phi_W4 = mid_way_f(Uno3List[j-1], u, dx, dt, G_C)

            G_DC = (Uno3List[0] - Uno3List[j])/dx
            G_CU = (Uno3List[j] - Uno3List[j-1])/dx

            if abs(G_DC - G_CU) <= 1.2 * abs(G_DC + G_CU)*0.5:
                G_C = G_DC - 2 * (dx + u * dt)/3*(G_DC - G_CU)/(2*dx)
            elif G_DC * G_CU > 0:
                G_C = sign(G_DC,0)*2*min(abs(G_DC), abs(G_CU))
            else:
                G_C = sign(G_DC,0)*min(abs(G_DC), abs(G_CU))

            phi_E4 = mid_way_f(Uno3List[j], u, dx, dt, G_C)

        else:
            #UNO2 plus

            phi_W = UnoPlusList[j - 1] + sign(UnoPlusList[j], UnoPlusList[j - 1]) * (dx - u * dt) * abs((UnoPlusList[j-1] - UnoPlusList[j - 2]) * (UnoPlusList[j] - UnoPlusList[j-1 ]) / (dx*dx))/( abs((UnoPlusList[j-1] - UnoPlusList[j-2]) / dx) + abs((UnoPlusList[j] - UnoPlusList[j-1]) / dx + 1e-308)) #* min(abs((temp_list[j] - temp_list[j - 1]) / dx), abs((temp_list[j + 1] - temp_list[j]) / dx))
            phi_E = UnoPlusList[j] + sign(UnoPlusList[j + 1], UnoPlusList[j]) * (dx - u * dt) * abs((UnoPlusList[j] - UnoPlusList[j - 1]) * (UnoPlusList[j + 1] - UnoPlusList[j]) / (dx*dx))/( abs((UnoPlusList[j] - UnoPlusList[j-1]) / dx) + abs((UnoPlusList[j+ 1] - UnoPlusList[j]) / dx + 1e-308))#min(abs((temp_list[j] - temp_list[j - 1]) / dx), abs((temp_list[j + 1] - temp_list[j]) / dx))

            #UNO2

            phi_W2 = UnoList[j - 1] + 0.5 * sign(UnoList[j], UnoList[j - 1]) * (dx - u * dt) * min(abs(UnoList[j-1] - UnoList[j - 2])/dx, abs(UnoList[j] - UnoList[j-1 ]) / dx) #* min(abs((temp_list[j] - temp_list[j - 1]) / dx), abs((temp_list[j + 1] - temp_list[j]) / dx))
            phi_E2 = UnoList[j] + 0.5 * sign(UnoList[j + 1], UnoList[j]) * (dx - u * dt) * min(abs((UnoList[j] - UnoList[j - 1]) / dx), abs((UnoList[j + 1] - UnoList[j]) / dx))
            
            #UNO3 minus

            G_DC = (Uno3MinusList[j] - Uno3MinusList[j-1])/dx
            G_CU = (Uno3MinusList[j-1] - Uno3MinusList[j-2])/dx

            if abs(G_DC - G_CU) <= 1.2 * abs(G_DC + G_CU)*0.5:
                G_C = G_DC - 2 * (dx + u * dt)/3*(G_DC - G_CU)/(2*dx)
            else:
                G_C = sign(G_DC,0)*2*abs(G_DC * G_CU)/(abs(G_DC) + abs(G_CU) + 1e-308)

            phi_W3 = mid_way_f(Uno3MinusList[j-1], u, dx, dt, G_C)

            G_DC = (Uno3MinusList[j+1] - Uno3MinusList[j])/dx
            G_CU = (Uno3MinusList[j] - Uno3MinusList[j-1])/dx

            if abs(G_DC - G_CU) <= 1.2 * abs(G_DC + G_CU)*0.5:
                G_C = G_DC - 2 * (dx + u * dt)/3*(G_DC - G_CU)/(2*dx)
            else:
                G_C = sign(G_DC,0)*2*abs(G_DC * G_CU)/(abs(G_DC) + abs(G_CU) + 1e-308)

            phi_E3 = mid_way_f(Uno3MinusList[j], u, dx, dt, G_C)

            #UNO3

            G_DC = (Uno3List[j] - Uno3List[j-1])/dx
            G_CU = (Uno3List[j-1] - Uno3List[j-2])/dx

            if abs(G_DC - G_CU) <= 1.2 * abs(G_DC + G_CU)*0.5:
                G_C = G_DC - 2 * (dx + u * dt)/3*(G_DC - G_CU)/(2*dx)
            elif G_DC * G_CU > 0:
                G_C = sign(G_DC,0)*2*min(abs(G_DC), abs(G_CU))
            else:
                G_C = sign(G_DC,0)*min(abs(G_DC), abs(G_CU))

            phi_W4 = mid_way_f(Uno3List[j-1], u, dx, dt, G_C)

            G_DC = (Uno3List[j+1] - Uno3List[j])/dx    
            G_CU = (Uno3List[j] - Uno3List[j-1])/dx

            if abs(G_DC - G_CU) <= 1.2 * abs(G_DC + G_CU)*0.5:
                
                G_C = G_DC - 2 * (dx + u * dt)/3*(G_DC - G_CU)/(2*dx)
            
            elif G_DC * G_CU > 0:
                
                G_C = sign(G_DC,0)*2*min(abs(G_DC), abs(G_CU))
            
            else:
                
                G_C = sign(G_DC,0)*min(abs(G_DC), abs(G_CU))

            phi_E4 = mid_way_f(Uno3List[j], u, dx, dt, G_C)


        # Update UnoPlusList using phi_W and phi_E
        temp_list[j] = temp_list[j] + (u * phi_W - u * phi_E) * dt / dx
        temp2_list[j] = temp2_list[j] + (u * phi_W2 - u * phi_E2) * dt / dx
        temp3_list[j] = temp3_list[j] + (u * phi_W3 - u * phi_E3) * dt / dx
        temp4_list[j] = temp4_list[j] + (u * phi_W4 - u * phi_E4) * dt / dx

    UnoPlusList = temp_list.copy()
    UnoList = temp2_list.copy()
    Uno3MinusList = temp3_list.copy()
    Uno3List = temp4_list.copy()

#print(UnoPlusList)
# Plot the list

rms = np.sqrt(np.mean((np.array(UnoPlusList) - np.array(square_list))**2))

print("RMS for UNO2+:", rms)

rms = np.sqrt(np.mean((np.array(UnoList) - np.array(square_list))**2))

print("RMS for UNO2:", rms)

rms = np.sqrt(np.mean((np.array(Uno3MinusList) - np.array(square_list))**2))

print("RMS for UNO3-:", rms)

rms = np.sqrt(np.mean((np.array(Uno3List) - np.array(square_list))**2))

print("RMS for UNO3:", rms)

plt.plot(UnoPlusList, label="UNO2+")
plt.plot(UnoList, label="UNO2")
plt.plot(Uno3MinusList, label="UNO3-")
plt.plot(Uno3List, label="UNO3")
plt.title("Square Function Plot")
plt.xlabel("Index")
plt.ylabel("Value")
plt.legend()
plt.show()
