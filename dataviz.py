import re
import numpy as np
import matplotlib.pyplot as plt

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