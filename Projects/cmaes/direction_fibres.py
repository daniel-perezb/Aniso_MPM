import cv2
import numpy as np
import matplotlib.pyplot as plt
from skimage.feature import structure_tensor, structure_tensor_eigenvalues
import sys
import math

# Load the image
image = cv2.imread('test1/maskedImage.png')
image_gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

# Load the mask
mask = cv2.imread('test1/initial_mask.png', cv2.IMREAD_GRAYSCALE)

# Normalize the image to [0, 1]
image_normalized = image_gray.astype(np.float32) / 255.0
mask = mask.astype(np.float32) / 255.0

# Structure Tensor Components
Axx, Axy, Ayy = structure_tensor(image_normalized, sigma=1)

# Build the full structure tensor
tensor = np.zeros(image_normalized.shape + (2, 2))
tensor[..., 0, 0] = Axx
tensor[..., 0, 1] = Axy
tensor[..., 1, 0] = Axy
tensor[..., 1, 1] = Ayy

# Compute eigenvalues
eigenvals = structure_tensor_eigenvalues(tensor)
lambda1, lambda2 = eigenvals[..., 0], eigenvals[..., 1]

# Compute Orientation
orientation = 0.5 * np.arctan2(2 * Axy, Axx - Ayy)
orientation_deg = (np.degrees(orientation) + 180) % 180  # Orientation in degrees [0, 180)

# Divide the Image into Cells
cell_size = 32  # Adjust cell size as needed
height, width = image_normalized.shape
rows = (height + cell_size - 1) // cell_size  # Ceiling division
cols = (width + cell_size - 1) // cell_size

# Prepare to store dominant orientations
cell_orientations = []

for i in range(rows):
    for j in range(cols):
        # Define the cell boundaries
        y_start = i * cell_size
        y_end = min((i + 1) * cell_size, height)
        x_start = j * cell_size
        x_end = min((j + 1) * cell_size, width)
        
        # Extract orientations within the cell
        cell_orientation = orientation_deg[y_start:y_end, x_start:x_end]
        
        # Flatten and remove NaN values (if any)
        angles = cell_orientation.flatten()
        angles = angles[~np.isnan(angles)]
        
        if len(angles) == 0:
            continue  # Skip empty cells
        
        # Convert angles to radians
        angles_rad = np.deg2rad(angles * 2)  # Multiply by 2 for orientations in [0, 180)
        
        # Compute mean direction
        sin_sum = np.sum(np.sin(angles_rad))
        cos_sum = np.sum(np.cos(angles_rad))
        mean_angle_rad = 0.5 * np.arctan2(sin_sum, cos_sum)
        mean_angle_deg = (np.rad2deg(mean_angle_rad) + 180) % 180  # Resulting angle in [0, 180)
        
        # Store the position and dominant orientation
        x_center = x_start + (x_end - x_start) / 2
        y_center = y_start + (y_end - y_start) / 2
        
        # Only store orientations if the center of the cell is inside the mask
        if mask[int(y_center), int(x_center)] > 0.5:  # Check mask value at the center
            cell_orientations.append({'position': (x_center, y_center), 'angle': mean_angle_deg})


# Define y_min and y_max for the filtered sub-image
y_min, y_max = 900, 1100 

fig, ax = plt.subplots(1, 2, figsize=(20, 10))

# Image in color with dominant orientations
ax[0].imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))
ax[0].axis('off')
ax[0].set_title('Dominant Fiber Orientations in Local Cells (Color Image)')

# Draw arrows representing the dominant orientation in each cell
arrow_length = cell_size / 2 

for item in cell_orientations:
    x_center, y_center = item['position']
    angle_deg = item['angle']
    theta_rad = np.deg2rad(angle_deg)
    dx = arrow_length * np.cos(theta_rad)
    dy = -arrow_length * np.sin(theta_rad)  # Negative due to image coordinates
    
    # Draw the arrow on the color image
    ax[0].arrow(x_center, y_center, dx, dy, color='red', width=1, head_width=5, length_includes_head=True)

# Second subplot: Filtered orientations within the y_min and y_max range
ax[1].imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))
ax[1].axis('off')
ax[1].set_title(f'Fiber Orientations Within y = [{y_min}, {y_max}]')

for item in cell_orientations:
    x_center, y_center = item['position']
    if y_min <= y_center <= y_max:
        angle_deg = item['angle']
        theta_rad = np.deg2rad(angle_deg)
        dx = arrow_length * np.cos(theta_rad)
        dy = -arrow_length * np.sin(theta_rad)
        
        # Draw the arrow only within the y range
        ax[1].arrow(x_center, y_center, dx, dy, color='blue', width=1, head_width=5, length_includes_head=True)

plt.tight_layout()
plt.show()

