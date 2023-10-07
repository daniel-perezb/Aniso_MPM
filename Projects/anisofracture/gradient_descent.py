import subprocess
import time
import os
import numpy as np

# Create a flag file to signal the C++ process to stop
flag_file_path = 'stop_flag.txt'
open(flag_file_path, 'w').close()

# Start the C++ code
cpp_process = subprocess.Popen(["./anisofracture", "-test", "8"])

# Define a function to check the condition
def check_condition():
    # Implement your condition-checking logic here
    # For example, you can run your point comparison or other code
    
    # Return True if condition is satisfied, False otherwise
    return True

def cleanup():
    # Remove the flag file to signal the C++ process to stop
    if os.path.exists(flag_file_path):
        os.remove(flag_file_path)

# Function to update values in 'parameters.txt'
def update_parameters_file():
    with open('parameters.txt', 'r') as file:
        lines = file.readlines()

    # Modify the values as needed based on your condition
    # For example, let's say you want to change Youngs to 25000 if the condition is not satisfied
    if not check_condition():
        for i in range(len(lines)):
            if 'Youngs:' in lines[i]:
                lines[i] = 'Youngs: 25000\n'

    # Write the updated contents back to the file
    with open('parameters.txt', 'w') as file:
        file.writelines(lines)


def video_variables():
    with open('/ziran2020/Data/TetMesh/mesh_files/visibles.txt', 'r') as file:
        visible_lines = file.readlines()
        values = [line.strip().split() for line in visible_lines]  # Use visible_lines

    # Convert values to float
    values = [[float(val) for val in row] for row in values]

    # Convert the list of lists to a NumPy array
    visible_array = np.array(values)  

    with open('/ziran2020/Data/TetMesh/mesh_files/data.txt', 'r') as file:
        data_lines = file.readlines()
        data_values = [line.strip().split() for line in data_lines]  # Use data_lines

    # Convert values to float
    data_values = [[float(val) for val in row] for row in data_values]

    # Convert the list of lists to a NumPy array
    data_array = np.array(data_values)  



# Open the data and variables
video_variables()

try:
    while True:
        if not check_condition():
            # Condition is not satisfied, so signal the C++ process to stop
            with open(flag_file_path, 'w') as flag_file:
                flag_file.write('stop')
            update_parameters_file()
            break

        # Optional: Add a delay to control how often you check the condition
        time.sleep(1)

except KeyboardInterrupt:
    # Handle Ctrl+C if needed
    pass
finally:
    cleanup()

# Wait for the C++ process to finish (if needed)
cpp_process.wait()
