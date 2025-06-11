import re
import numpy as np
import matplotlib.pyplot as plt

# Paste your data here as a multiline string
data_str = """
6.66192e-078
-2.36844e-077
-2.16491e-076
 6.81436e-076
 2.53773e-075
 2.07196e-074
 4.93986e-074
 4.41541e-073
-1.60712e-072
-1.53307e-071
-2.48233e-071
-3.54733e-071
-6.11902e-070
-4.76656e-070
  9.9985e-069
  6.9092e-068
 1.97043e-067
 7.58198e-067
 2.54257e-066
 7.23049e-066
-1.28155e-065
-2.07292e-064
-1.31218e-063
-1.05664e-062
-2.70463e-062
 -1.5317e-062
 4.71369e-061
 2.08667e-060
 6.84542e-060
 4.95802e-059
 1.35016e-058
 4.65878e-058
 1.28918e-057
 3.48252e-057
 8.73898e-057
 6.38865e-057
 2.50509e-055
 5.92086e-055
 2.39121e-054
-6.42197e-054
 4.49873e-053
 5.34355e-052
-1.87689e-051
-1.80275e-050
 3.05064e-052
 2.08422e-051
 9.69947e-050
-2.14604e-048
 -8.6826e-048
-5.33001e-046
 -1.7778e-045
 1.02786e-044
-7.78216e-044
-2.91563e-043
-9.09685e-043
 5.52924e-042
 1.83694e-041
 1.16159e-040
 1.87497e-040
 2.23554e-039
-2.74411e-041
 5.25386e-038
 2.57628e-037
 9.99288e-037
 3.51489e-036
 1.17735e-035
 3.83168e-035
 1.22563e-034
 3.87644e-034
 1.21704e-033
 3.80152e-033
 1.18319e-032
 3.67274e-032
 1.13772e-031
 3.51845e-031
 1.08654e-030
 3.35109e-030
 1.03231e-029
 3.17643e-029
 9.76312e-029
 2.99753e-028
 9.19318e-028
 2.81639e-027
 8.61868e-027
 2.63454e-026
 8.04414e-026
 2.45334e-025
 7.47367e-025
 2.27405e-024
 6.91111e-024
 2.09783e-023
 6.36004e-023
 1.92579e-022
 5.82379e-022
 1.75892e-021
 5.30539e-021
 1.59813e-020
 4.80753e-020
 1.44423e-019
 4.33256e-019
 1.29789e-018
 3.88245e-018
 1.15968e-017
 3.45876e-017
 1.03002e-016
 3.06268e-016
  9.0923e-016
 2.69496e-015
 7.97488e-015
   2.356e-014
 6.94855e-014
 2.04581e-013
 6.01278e-013
 1.76404e-012
 5.16597e-012
 1.51004e-011
 4.40558e-011
 1.28286e-010
 3.72817e-010
 1.08129e-009
 3.12963e-009
 9.03928e-009
 2.60522e-008
  7.4921e-008
 2.14976e-007
 6.15437e-007
 1.75777e-006
 5.00839e-006
 1.42355e-005
 4.03605e-005
  0.000114137
  0.000321923
  0.000905542
   0.00254018
   0.00710542
    0.0198177
    0.0551085
     0.152776
     0.422206
      1.16302
      3.19307
      8.73658
      23.8202
      64.7107
      175.139
      472.191
      1268.02
      3391.21
      9031.14
      23945.8
        63205
       166050
       434128
 1.12931e+006
  2.9224e+006
 7.52169e+006
 1.92506e+007
 4.89812e+007
  1.2387e+008
 3.11271e+008
 7.77016e+008
 1.92624e+009
 4.74065e+009
 1.15789e+010
 2.80564e+010
 6.74152e+010
 1.60567e+011
 3.78893e+011
 8.85349e+011
 2.04739e+012
 4.68274e+012
 1.05855e+013
 2.36321e+013
 5.20592e+013
 1.13052e+014
 2.41756e+014
 5.08466e+014
 1.05033e+015
 2.12756e+015
 4.21823e+015
 8.16887e+015
 1.54143e+016
 2.82629e+016
 5.01974e+016
 8.60649e+016
 1.41928e+017
 2.24304e+017
 3.38634e+017
 4.87228e+017
  6.6753e+017
 8.71921e+017
 1.08995e+018
 1.31279e+018
 1.53879e+018
  1.7795e+018
 2.06903e+018
 2.49205e+018
 3.31508e+018
 6.45262e+018
 9.66082e+019
 9.99824e+019
 9.99969e+019
 9.99962e+019
 9.99953e+019
 9.99945e+019
 9.99937e+019
  9.9993e+019
 9.99923e+019
 9.99917e+019
 9.99911e+019
 9.99907e+019
 9.99902e+019
 9.99899e+019
 9.99897e+019
 9.99895e+019
 9.99894e+019
 9.99894e+019
 9.99895e+019
 9.99896e+019
 9.99899e+019
 9.99901e+019
 9.99905e+019
 9.99909e+019
 9.99914e+019
 9.99919e+019
 9.99925e+019
 9.99932e+019
 9.99939e+019
 9.99947e+019
 9.99955e+019
 9.99964e+019
 9.99973e+019
 9.99983e+019
 9.99993e+019
       1e+020
 1.00001e+020
 1.00003e+020
 1.00004e+020
 1.00005e+020
 1.00006e+020
 1.00007e+020
 1.00008e+020
  1.0001e+020
 1.00011e+020
 1.00012e+020
 1.00013e+020
 1.00014e+020
 1.00015e+020
 1.00016e+020
 1.00016e+020
 1.00017e+020
 1.00017e+020
 1.00018e+020
 1.00018e+020
 1.00018e+020
 1.00018e+020
 1.00018e+020
 1.00018e+020
 1.00018e+020
 1.00017e+020
 1.00017e+020
 1.00016e+020
 1.00015e+020
 1.00015e+020
 1.00014e+020
 1.00013e+020
 1.00012e+020
 1.00011e+020
  1.0001e+020
 1.00009e+020
 1.00007e+020
 1.00006e+020
 1.00005e+020
 1.00004e+020
 1.00003e+020
 1.00002e+020
       1e+020
 9.99993e+019
 9.99983e+019
 9.99973e+019
 9.99963e+019
 9.99954e+019
 9.99946e+019
 9.99939e+019
 9.99933e+019
 9.99928e+019
 9.99925e+019
 9.99923e+019
 9.99922e+019
 9.99923e+019
 9.99925e+019
 9.99929e+019
 9.99934e+019
  9.9994e+019
 9.99948e+019
 9.99958e+019
 9.99967e+019
  9.9981e+019
 9.62008e+019
 6.47516e+018
 3.35142e+018
  2.4017e+018
  1.7207e+018
 1.13149e+018
 6.60945e+017
 3.40187e+017
 1.55154e+017
 6.34756e+016
 2.36205e+016
 8.09883e+015
 2.58682e+015
 7.76637e+014
 2.20767e+014
  5.9772e+013
 1.54896e+013
 3.85786e+012
   9.267e+011
 2.15337e+011
 4.85307e+010
 1.06321e+010
 2.26881e+009
 4.72419e+008
 9.61378e+007
 1.91477e+007
 3.73728e+006
       715670
       134601
      24887.4
      4527.83
      811.208
      143.229
      24.9392
      4.28514
     0.727005
     0.121853
    0.0201876
   0.00330744
  0.000536106
 8.60089e-005
 1.36629e-005
 2.14984e-006
 3.35187e-007
 5.17998e-008
 7.93713e-009
  1.2062e-009
 1.81853e-010
 2.72066e-011
 4.04012e-012
 5.95636e-013
 8.72031e-014
 1.26806e-014
 1.83187e-015
 2.62954e-016
 3.75124e-017
 5.31934e-018
 7.49897e-019
 1.05118e-019
 1.46538e-020
 2.03182e-021
  2.8025e-022
 3.84583e-023
 5.25138e-024
 7.13596e-025
 9.65113e-026
 1.29928e-026
 1.74129e-027
 2.32345e-028
 3.08697e-029
 4.08425e-030
 5.38164e-031
 7.06285e-032
 9.23309e-033
 1.20241e-033
 1.56004e-034
 2.01664e-035
 2.59755e-036
 3.33409e-037
 4.26481e-038
 5.43702e-039
 6.90863e-040
 8.75026e-041
 1.10478e-041
 1.39054e-042
  1.7449e-043
 2.18304e-044
 2.72322e-045
 3.38732e-046
 4.20152e-047
 5.19704e-048
 6.41101e-049
 7.88751e-050
  9.6787e-051
 1.18461e-051
 1.44624e-052
 1.76127e-053
  2.1397e-054
 2.59322e-055
 3.13547e-056
 3.78233e-057
 4.55226e-058
 5.46666e-059
 6.55028e-060
  7.8317e-061
 9.34387e-062
 1.11246e-062
 1.32174e-063
 1.56721e-064
 1.85454e-065
 2.19024e-066
 2.58168e-067
 3.03728e-068
 3.56656e-069
 4.18032e-070
 4.89076e-071
 5.71165e-072
 6.65851e-073
 7.74878e-074
 9.00207e-075
 1.04403e-075
 1.20881e-076
  1.3973e-077
 1.61254e-078
 1.85796e-079
 2.13735e-080
 2.45492e-081
 2.81534e-082
 3.22377e-083
 3.68594e-084
 4.20814e-085
 4.79732e-086
 5.46112e-087
 6.20796e-088
 7.04705e-089
 7.98849e-090
 9.04335e-091
 1.02237e-091
 1.15427e-092
 1.30148e-093
 1.46556e-094
 1.64821e-095
 1.85128e-096
 2.07677e-097
 2.32685e-098
 2.60386e-099
 2.91033e-100
 3.24901e-101
 3.62284e-102
 4.03498e-103
 4.48884e-104
 4.98808e-105
 5.53663e-106
 6.13869e-107
 6.79875e-108
 7.52163e-109
 8.31245e-110
 9.17669e-111
 1.01202e-111
 1.11491e-112
   1.227e-113
   1.349e-114
 1.48164e-115
 1.62571e-116
 1.78205e-117
 1.95152e-118
 2.13506e-119
 2.33364e-120
 2.54831e-121
 2.78013e-122
 3.03026e-123
 3.29989e-124
 3.59028e-125
 3.90275e-126
 4.23867e-127
  4.5995e-128
 4.98673e-129
 5.40195e-130
 5.84679e-131
 6.32297e-132
 6.83227e-133
 7.37654e-134
 7.95769e-135
 8.57773e-136
 9.23872e-137
 9.94281e-138
 1.06922e-138
 1.14892e-139
 1.23361e-140
 1.32355e-141
 1.41897e-142
 1.52015e-143
 1.62734e-144
 1.74081e-145
 1.86086e-146
 1.98775e-147
  2.1218e-148
  2.2633e-149
 2.41254e-150
  2.5713e-151
 2.33185e-152
"""  # Replace "..." with the full data list

# Convert the string data into a NumPy array
data = np.array([float(line) for line in data_str.strip().splitlines()])

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(data)
#plt.yscale('log')  # Optional: log scale to handle wide dynamic range
plt.title("Data Visualization")
plt.xlabel("Index")
plt.ylabel("Value")
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.show()


# Define times to look for (in seconds)
target_times = [0, 1.09067e-007, 4.95446e-007, 1.01406e-006, 2.25037e-006, 5.00491e-006, 7.50093e-006]
tolerance = 1e-7

filename = '../output 10 micro sec 9000s/e_dens.txt'


# Read file content
with open(filename, 'r') as file:
    content = file.read()

# Split by blocks
blocks = content.split('----')

times = []
data_blocks = []

for block in blocks:
    lines = block.strip().split('\n')
    #print(block)
    if len(lines) < 2:
        continue  # Skip if not enough lines

    # Extract time
    time_match = re.search(r'Time:\s*([0-9eE\+\-\.]+)', lines[0])
    print(time_match)
    if not time_match:
        continue

    time = float(time_match.group(1))

    # Parse data line (single line of numbers)
    try:
        data_line = lines[1]
        data = list(map(float, data_line.split()))
        times.append(time)
        data_blocks.append(data)
    except Exception as e:
        print(f"Error parsing data at time {time:.2e}: {e}")

# Convert to numpy arrays
times = np.array(times)
data_blocks = np.array(data_blocks, dtype=object)  # allow ragged arrays just in case

# Generate 102 values from 0 to 0.0003
x = np.linspace(0, 0.0003, 102)

# Plot values at requested times
for t in target_times:
    idx = np.where(np.abs(times - t) < tolerance)[0]
    if idx.size > 0:
        i = idx[0]
        plt.plot(x,data_blocks[i], label=f'Time = {times[i]:.2e} s')
    else:
        print(f"Time {t:.2e} s not found within tolerance.")

plt.xlabel("Z (m)")
plt.ylabel("m^-3")
plt.legend()
plt.title("electron density ")
plt.grid(True)
plt.show()

# File containing flux data
filename = '../output 10 micro sec 9000s/current_dens.txt'  # replace with your file path

# Initialize lists to store flux values
left_flux = []
right_flux = []

# Read the file
with open(filename, 'r') as file:
    for line in file:
        # Strip any extra spaces or newline characters
        line = line.strip()
        
        # Split the line into left and right flux values (separated by '|')
        flux_values = line.split('|')
        
        # Convert the flux values from string to float
        left_flux.append(float(flux_values[0].strip()))  # Left flux (particles/s)
        right_flux.append(float(flux_values[1].strip()))  # Right flux (particles/s)

# Convert flux data into numpy arrays for easier manipulation
left_flux = np.array(left_flux)
right_flux = np.array(right_flux)

left_flux = left_flux[1:]
right_flux = right_flux[1:]

# Compute net particle flux and current
net_flux = left_flux + right_flux  # sum left + right fluxes for each time step
current = net_flux  # Compute the current in amperes (A)

# Create time steps based on the number of flux values (you can adjust as needed)
time_steps = np.arange(1, len(left_flux) + 1)

# Plot all on the same graph
plt.figure(figsize=(10, 6))

# Plot Left flux
plt.plot(times[2:], left_flux, label='Left End Flux', color='blue')

# Plot Right flux
plt.plot(times[2:], right_flux, label='Right End Flux', color='red')

# Plot Current
plt.plot(times[2:], current, label='Current', color='purple')

# Add labels and title
plt.xlabel('Time (s))')
plt.ylabel('Current (A)')
plt.title('Current')

# Show the legend
plt.legend()

# Show the plot
plt.show()

# File containing potential and dt data
filename = '../output 10 micro sec 9000s/time_steps.txt'  # replace with your file path

# Initialize lists to store dt values
dt_values = []

# Read the file and extract dt values
with open(filename, 'r') as file:
    for line in file:
        # Strip any extra spaces or newline characters
        line = line.strip()
        
        # Extract the 'dt' value from the line using string splitting
        if 'dt' in line:
            # Split the line at 'dt' and take the second part, then extract the numerical value
            dt_value = line.split('dt:')[1].strip()
            dt_values.append(float(dt_value))

# Convert dt values into numpy array for easier manipulation
dt_values = np.array(dt_values)

# Create iteration numbers (1-based index for each time step)
iterations = np.arange(1, len(dt_values) + 1)

# Plot dt across iterations
plt.figure(figsize=(10, 6))
plt.plot(times[1:], dt_values, color='blue')

# Add labels and title
plt.xlabel('Iteration')
plt.ylabel('dt (seconds)')
plt.title('Time Step (dt)')

# Show the legend
plt.legend()

# Show the plot
plt.show()