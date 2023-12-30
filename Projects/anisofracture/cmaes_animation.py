import subprocess
import os
import time
import glob
import csv
import pickle
import time
import numpy as np
from cmaes import CMA

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

    def run(self, timeout_seconds=60):
        self.clean_output_directory()
        self.start_cpp_process()
        print("Started C++ process")
        
        try:
            self.cpp_process.wait(timeout=timeout_seconds * 60 * 8)  # Convert minutes to seconds
            print("C++ process finished")
            
        except subprocess.TimeoutExpired:
            print("C++ process timed out after 8 hours")
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
        # Loop through all files in the output folder
        for file_name in os.listdir(self.output_folder):
            file_path = os.path.join(self.output_folder, file_name)
            os.remove(file_path)
  

def read_parameters(filepath="/Fibre_directions/Aniso_MPM/Projects/anisofracture/parameters.txt"):
    params = {}
    with open(filepath, 'r') as f:
        for line in f:
            key, value = line.split(":")
            params[key.strip()] = float(value.strip())
    # Return values in the order of initial mean
    return [params["Youngs"], params["nu"], params["rho"], params["eta"], params["percentage"], params["fiber"]]

def write_parameters(params, filepath="/Fibre_directions/Aniso_MPM/Projects/anisofracture/parameters.txt"):
    with open(filepath, 'w') as f:
        for key, value in params.items():
            f.write("{}: {}\n".format(key, value))

def normalize(params, lower_bounds, upper_bounds):
    return (params - lower_bounds) / (upper_bounds - lower_bounds)

def denormalize(normalized_params, lower_bounds, upper_bounds):
    return normalized_params * (upper_bounds - lower_bounds) + lower_bounds

def compute_mse():
    subprocess.run(["calculate_mse/build/./open"])
    with open('mse.txt', 'r') as file:
        mse = file.read()

    mse = float(mse)
    return mse

def read_parameters_from_csv(loop):
    with open('solutions_mse.csv', 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header row
        for index, row in enumerate(reader):
            if index == loop:
                # Splitting the string by spaces instead of commas
                solution_str = row[0].strip('[]').split()
                solution = np.array([float(x.strip()) for x in solution_str])
                return solution
    return None


def objective_function(normalized_params_array, loop):
    params_from_csv = read_parameters_from_csv(loop)

    if params_from_csv is not None:
        # Convert the parameters back to their original scale
        params_array = denormalize(params_from_csv, lower_bounds, upper_bounds)
        params_dict = {"Youngs": params_array[0], "nu": params_array[1], "rho": params_array[2], "eta": params_array[3], "percentage": params_array[4], "fiber": params_array[5]}
        write_parameters(params_dict)
        subprocess.run(["cp", "-R", f"/home/daniel/Fleece_animations/Extract_STL/folder_{loop+1}/RIG_fleece_fibres/data_0.dat", "output/RIG_fleece_fibres/"])
        subprocess.run(["cp", "-R", f"/home/daniel/Fleece_animations/Extract_STL/folder_{loop+1}/RIG_fleece_fibres/data_70.dat", "output/RIG_fleece_fibres/"])
        mse = compute_mse()
        print(f"Computed MSE: {mse}")
        time.sleep(3)
    
    else:
        # If not, proceed with running the animation and computing MSE
        params_array = denormalize(normalized_params_array, lower_bounds, upper_bounds)
        params_dict = {"Youngs": params_array[0], "nu": params_array[1], "rho": params_array[2], "eta": params_array[3], "percentage": params_array[4], "fiber": params_array[5]}
        write_parameters(params_dict)
        success = animation.run()

        if not success:
            return float('inf')  

        mse = compute_mse()
    
    return mse


if __name__ == '__main__':
    
    animation = Animation()

    initial_mean = np.array(read_parameters())
    lower_bounds = np.array([100, 0.25, 1, 0.01, .1, 10])
    upper_bounds = np.array([500000, .45, 800, 0.5, 9999, 300])
    
    normalized_initial_mean = normalize(initial_mean, lower_bounds, upper_bounds)
    bounds_array = np.vstack([np.zeros_like(lower_bounds), np.ones_like(upper_bounds)]).T
    optimizer = CMA(mean=normalized_initial_mean, sigma=1.3, bounds=bounds_array)
    
    print("Created a new optimizer")

    all_solutions = []
    all_mse_values = []
    loop = 0

    for generation in range(5):
        solutions = []

        for _ in range(optimizer.population_size):
            solution = optimizer.ask()
            mse = objective_function(solution, loop)
            solutions.append((solution, mse))
            
            # Save the solution and its mse
            all_solutions.append(solution)
            all_mse_values.append(mse)
            
            loop += 1
            
            folder_name = f"/home/daniel/Fleece_animations/Extract_STL/folder_{loop}"

            if not os.path.exists(folder_name):
                os.makedirs(folder_name, exist_ok=True)
                # Copy the contents from the source to the newly created folder 
                subprocess.run(["cp", "-R", "/Fibre_directions/Aniso_MPM/Projects/anisofracture/output/RIG_fleece_fibres/", folder_name])

            # Save the optimizer's state
            with open('/Fibre_directions/Aniso_MPM/Projects/anisofracture/optimizer_state.pkl', 'wb') as f:
                pickle.dump(optimizer, f)

            # Save all solutions and their MSE values to a CSV file
            with open('new_solutions_mse.csv', 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(['Solution', 'MSE'])
                for solution, mse in zip(all_solutions, all_mse_values):
                    writer.writerow([solution, mse])

        
        optimizer.tell(solutions)


    print("All solutions and their MSE values have been saved to solutions_mse.csv")