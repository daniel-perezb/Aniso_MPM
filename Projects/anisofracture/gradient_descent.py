import subprocess
import os
import time
import numpy as np
from mse import calculate_mse

class Animation:
    def __init__(self):
        self.cpp_process = None
        self.output_folder = "/Fibre_directions/Aniso_MPM/Projects/anisofracture/output/RIG_fleece_fibres/"
        self.target_file = "data_{}.dat"  # Specify the target file pattern
        self.target_file_number = 0  # Initialize the target file number

    def start_cpp_process(self):
        self.target_file = "data_{}.dat".format(self.target_file_number)
        self.cpp_process = subprocess.Popen([
            "/Fibre_directions/Aniso_MPM/Projects/anisofracture/./anisofracture",
            "-test", "8"
        ])
        # self.target_file_number += 2  # Increment the target file number

    def stop_cpp_process(self):
        if self.cpp_process and self.cpp_process.poll() is None:
            self.cpp_process.terminate()
    
    def run(self):
        self.start_cpp_process()
        print("Started C++ process")
        self.cpp_process.wait()
        print("C++ process finished")

'''
    def run(self):
        run_code = True
        while run_code:
            self.start_cpp_process()
            print("Started C++ process for {}".format(self.target_file))
            while True:
                # Check if the target file exists
                if os.path.exists(os.path.join(self.output_folder, self.target_file)):
                    message = 'Target file "{}" found. Stopping the C++ process...'.format(self.target_file)
                    print(message)
                    self.stop_cpp_process()
                    break  # Exit the inner loop when the target file is found
                    run_code = False
                time.sleep(0.5)  # Adjust the polling interval as needed
'''

def read_parameters(filepath="/Fibre_directions/Aniso_MPM/Projects/anisofracture/parameters.txt"):
    params = {}
    with open(filepath, 'r') as f:
        for line in f:
            key, value = line.split(":")
            params[key.strip()] = float(value.strip())
    return params

def write_parameters(params, filepath="/Fibre_directions/Aniso_MPM/Projects/anisofracture/parameters.txt"):
    with open(filepath, 'w') as f:
        for key, value in params.items():
            f.write("{}: {}\n".format(key, value))

def compute_mse():
    dataFilePath = "/Fibre_directions/Aniso_MPM/Data/TetMesh/mesh_files/data.txt"
    logicFilePath = "/Fibre_directions/Aniso_MPM/Data/TetMesh/mesh_files/visibles.txt"
    scaleFilePath = "/Fibre_directions/Aniso_MPM/Data/TetMesh/mesh_files/scale_values.txt"
    directory = "/Fibre_directions/Aniso_MPM/Projects/anisofracture/output/RIG_fleece_fibres/"
    mse = calculate_mse(dataFilePath, logicFilePath, scaleFilePath, directory)
    return mse

def compute_gradient(animation, epsilon=1e-5):
    initial_params = read_parameters()
    initial_mse = compute_mse()
    
    gradients = {}
    for param in initial_params:
        perturbed_params = initial_params.copy()
        perturbed_params[param] += epsilon
        write_parameters(perturbed_params)
        animation.run()
        perturbed_mse = compute_mse()

        gradient = (perturbed_mse - initial_mse) / epsilon
        gradients[param] = gradient

    return gradients

if __name__ == '__main__':
    animation = Animation()
    learning_rate = .001
    iterations = 10
    mse_history = []


    # Define scaling factors for each parameter
    
    scaling_factors = {
        'Youngs': 50000000,
        'nu': 500,
        'rho': 200000,
        'residual_stress': 10
    }

    for _ in range(iterations):
        gradients = compute_gradient(animation)
        params = read_parameters()
        for param, gradient in gradients.items():
            params[param] -= learning_rate * gradient
            #update_amount = learning_rate * gradient * scaling_factors[param]
            #params[param] -= update_amount
        write_parameters(params)
        mse = compute_mse()
        mse_history.append(mse)
        with open('/Fibre_directions/Aniso_MPM/Projects/anisofracture/all_mse.txt', 'w') as f:
        for grad_value in mse_history:
            f.write(str(mse_value) + '\n')
    

    
    with open('/Fibre_directions/Aniso_MPM/Projects/anisofracture/all_mse.txt', 'w') as f:
        for mse_value in mse_history:
            f.write(str(mse_value) + '\n')
    
    print("Yet")
    
    '''
    animation.run()
    print("Current MSE:", mse)
    print("-----")
    '''
    