import subprocess
import os
import time
import glob
import csv
import pickle
import time
import numpy as np
from cmaes import CMA
from image_comparison import *

# Use pixel similarity for cmaes
image_error = False

class Animation:
    def __init__(self):
        self.cpp_process = None
        self.output_folder = "../anisofracture/output/RIG_fleece_fibres/"
        self.target_file_pattern = "data_{}.dat"  # Specify the target file pattern
        self.target_file_number = 2  # Initialize the target file number

    def start_cpp_process(self):
        self.target_file = self.target_file_pattern.format(self.target_file_number)
        self.cpp_process = subprocess.Popen([
            "../anisofracture/./anisofracture",
            "-test", "27"
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
  

def read_parameters(filepath="parameters.txt"):
    params = {}
    with open(filepath, 'r') as f:
        for line in f:
            if ':' in line: 
                print(line)
                key, value = line.split(":")
                params[key.strip()] = float(value.strip())
    # Return values in the order of initial mean
    return [params["Youngs"], params["nu"], params["rho"], params["eta"], params["percentage"], params["fiber"]]

def write_parameters(params, filepath="parameters.txt"):
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

def objective_function(normalized_params_array, loop):

    # Convert the parameters back to their original scale
    params_array = denormalize(normalized_params_array, lower_bounds, upper_bounds)
    params_dict = {"Youngs": params_array[0], "nu": params_array[1], "rho": params_array[2], "eta": params_array[3], "percentage": params_array[4], "fiber": params_array[5]}
    write_parameters(params_dict)
    
    success = animation.run()

    mse = compute_mse()
    print(f"Computed MSE: {mse}")
    time.sleep(3)
        
    if not success:
        return float('inf')  

    mse = compute_mse()
    
    return mse

def compute_pixel_error(loop): 
   
    target_image_path = 'skirted_fleece.png'
    scale = 2.1
    translation = (-160, 90)
    rect_coords = (350, 600, 220, 220)  # (x, y, width, height)
    screenshot_coords = (1470, 450, 1920, 1100)  # (x1, y1, x2, y2)
    percentage_diff = compute_percentage_difference(target_image_path, scale, translation, rect_coords, screenshot_coords)
    print(f"Computed Pixel Error: {percentage_diff}")
    return percentage_diff


if __name__ == '__main__':
    
    animation = Animation()

    initial_mean = np.array(read_parameters())
    lower_bounds = np.array([100, 0.25, 1, 0.01, .1, 10])
    upper_bounds = np.array([500000, .45, 800, 0.5, 9999, 300])
    
    normalized_initial_mean = normalize(initial_mean, lower_bounds, upper_bounds)
    bounds_array = np.vstack([np.zeros_like(lower_bounds), np.ones_like(upper_bounds)]).T
    
    if os.path.exists('optimizer_state.pkl'):
        with open('optimizer_state.pkl', 'rb') as f:
            optimizer = pickle.load(f)
        print("Loaded existing optimizer state")
    else:
        optimizer = CMA(mean=normalized_initial_mean, sigma=1.3, bounds=bounds_array)
        print("Created a new optimizer")

    print("Created a new optimizer")

    all_solutions = []
    all_mse_values = []
    all_pixel_error = []
    loop = 0

    for generation in range(10):
        solutions = []

        for _ in range(optimizer.population_size):
            solution = optimizer.ask()
            mse = objective_function(solution, loop)
            if mse == float('inf'):
                pixel_error = float('inf')
                print(f"Computed Pixel Error: {pixel_error}")
                error: float('inf')
            else:
                if image_error:
                    pixel_error = compute_pixel_error(loop)
                    error =  (.5 * pixel_error) + (400 * mse)
                else:
                    error = 800 * mse

            solutions.append((solution, error))
            
            # Save the solution and its mse
            all_solutions.append(solution)
            all_mse_values.append(mse)
            all_pixel_error.append(pixel_error)
            
            loop += 1
            
            folder_name = f"../../Data/Simulations/folder_{loop}"

            if not os.path.exists(folder_name):
                os.makedirs(folder_name, exist_ok=True)
                # Copy the contents from the source to the newly created folder 
                # subprocess.run(["cp", "-R", "output/RIG_fleece_fibres/", folder_name])
                # subprocess.run(["cp", "-R", "output/RIG_fleece_fibres/partio_7.bgeo", "../../Data/final/"])
        
            # Save the optimizer's state
            with open('optimizer_state.pkl', 'wb') as f:
                pickle.dump(optimizer, f)

            # Save all solutions and their MSE values to a CSV file
            with open('solutions_mse.csv', 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(['Solution', 'MSE','Pixel'])
                for solution, mse, pixel_error in zip(all_solutions, all_mse_values, all_pixel_error):
                    writer.writerow([solution, mse, pixel_error])

        optimizer.tell(solutions)


    print("All solutions and their MSE values have been saved to solutions_mse.csv")