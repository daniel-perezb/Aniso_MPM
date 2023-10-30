import subprocess
import os
import time
import glob
import csv
import pickle
import time
import numpy as np
from cmaes import CMA
from mse import calculate_mse

class Animation:
    def __init__(self):
        self.cpp_process = None
        self.output_folder = "/Fibre_directions/Aniso_MPM/Projects/anisofracture/output/RIG_fleece_fibres/"
        self.target_file_pattern = "data_{}.dat"  # Specify the target file pattern
        self.target_file_number = 2  # Initialize the target file number

    def start_cpp_process(self):
        self.target_file = self.target_file_pattern.format(self.target_file_number)
        self.cpp_process = subprocess.Popen([
            "/Fibre_directions/Aniso_MPM/Projects/anisofracture/./anisofracture",
            "-test", "8"
        ])

    def stop_cpp_process(self):
        if self.cpp_process and self.cpp_process.poll() is None:
            self.cpp_process.terminate()

    def run(self, timeout_minutes=60):
        self.clean_output_directory()
        self.start_cpp_process()
        print("Started C++ process")
        
        try:
            self.cpp_process.wait(timeout=timeout_minutes * 60 * 4)  # Convert minutes to seconds
            print("C++ process finished")
            
        except subprocess.TimeoutExpired:
            print("C++ process timed out after {} minutes".format(timeout_minutes))
            self.stop_cpp_process()
            return False
        
        if not self.output_file_exists():
            print("Output file not found: {}".format(self.target_file))
            return False
        return True


    def output_file_exists(self):
        output_file_path = os.path.join(self.output_folder, self.target_file)
        return os.path.exists(output_file_path)
    
    def clean_output_directory(self):
        files_to_remove = glob.glob(os.path.join(self.output_folder, "data_*.dat"))
        for file_path in files_to_remove:
            os.remove(file_path)
    

def read_parameters(filepath="/Fibre_directions/Aniso_MPM/Projects/anisofracture/parameters.txt"):
    params = {}
    with open(filepath, 'r') as f:
        for line in f:
            key, value = line.split(":")
            params[key.strip()] = float(value.strip())
    # Return values in the order of initial mean
    return [params["Youngs"], params["nu"], params["rho"], params["residual_stress"], params["eta"], params["cfl"], params["flip_pic_ratio"], params["percentage"]]

def write_parameters(params, filepath="/Fibre_directions/Aniso_MPM/Projects/anisofracture/parameters.txt"):
    with open(filepath, 'w') as f:
        for key, value in params.items():
            f.write("{}: {}\n".format(key, value))

def normalize(params, lower_bounds, upper_bounds):
    return (params - lower_bounds) / (upper_bounds - lower_bounds)

def denormalize(normalized_params, lower_bounds, upper_bounds):
    return normalized_params * (upper_bounds - lower_bounds) + lower_bounds

def compute_mse():
    dataFilePath = "/Fibre_directions/Aniso_MPM/Data/TetMesh/mesh_files/data.txt"
    logicFilePath = "/Fibre_directions/Aniso_MPM/Data/TetMesh/mesh_files/visibles.txt"
    scaleFilePath = "/Fibre_directions/Aniso_MPM/Data/TetMesh/mesh_files/scale_values.txt"
    directory = "/Fibre_directions/Aniso_MPM/Projects/anisofracture/output/RIG_fleece_fibres/"
    mse = calculate_mse(dataFilePath, logicFilePath, scaleFilePath, directory)
    return mse

def objective_function(normalized_params_array):
    params_array = denormalize(normalized_params_array, lower_bounds, upper_bounds)
    params_dict = {"Youngs": params_array[0], "nu": params_array[1], "rho": params_array[2], "residual_stress": params_array[3], "eta": params_array[4], "cfl": params_array[5], "flip_pic_ratio": params_array[6], "percentage": params_array[7]}
    write_parameters(params_dict)
    success = animation.run()

    # Return a high error value to indicate failure if animation does not run
    if not success:
        return float('inf')  
 
    mse = compute_mse()
    return mse


if __name__ == '__main__':
    animation = Animation()

    initial_mean = np.array(read_parameters())
    lower_bounds = np.array([200, 0.25, 1, 0.0001, 0.01, 0.01, 0, .2])
    upper_bounds = np.array([5000000, .45, 800, 0.9, 0.5, 0.5, 10, 99999])

    normalized_initial_mean = normalize(initial_mean, lower_bounds, upper_bounds)
    bounds_array = np.vstack([np.zeros_like(lower_bounds), np.ones_like(upper_bounds)]).T
    optimizer = CMA(mean=normalized_initial_mean, sigma=1.3, bounds=bounds_array)
    print("Created a new optimizer")

    all_solutions = []
    all_mse_values = []

    for generation in range(100):
        solutions = []
        for _ in range(optimizer.population_size):
            solution = optimizer.ask()
            mse = objective_function(solution)
            solutions.append((solution, mse))
            
            # Save the solution and its mse
            all_solutions.append(solution)
            all_mse_values.append(mse)

        optimizer.tell(solutions)

        # Save the optimizer's state
        with open('/Fibre_directions/Aniso_MPM/Projects/anisofracture/optimizer_state.pkl', 'wb') as f:
            pickle.dump(optimizer, f)

    # Save all solutions and their MSE values to a CSV file
    with open('solutions_mse.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Solution', 'MSE'])
        for solution, mse in zip(all_solutions, all_mse_values):
            writer.writerow([solution, mse])

    print("All solutions and their MSE values have been saved to solutions_mse.csv")