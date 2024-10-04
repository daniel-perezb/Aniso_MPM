import cv2
import numpy as np
import matplotlib.pyplot as plt


def load_image(image_path):
    """Load an image in grayscale"""
    return cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)

def fill_holes(image):
  """Fill small holes in the given binary image based on a maximum hole size threshold."""
  # Invert the image to treat holes as foreground
  inverted_image = cv2.bitwise_not(image)

  # Find connected components (holes) in the inverted image
  num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(inverted_image, connectivity=4)

  # Create an empty mask to fill holes based on their size
  filled_image = np.zeros_like(image)
  
  # Define a max threshold for hold size
  max_hole_size=80

  # Iterate through the connected components
  for i in range(1, num_labels):
      area = stats[i, cv2.CC_STAT_AREA]
        
      # Fill only if the hole size is less than or equal to the max_hole_size
      if area <= max_hole_size:
          filled_image[labels == i] = 255

  # Combine the filled holes with the original image
  result_image = cv2.bitwise_or(image, filled_image)

  return result_image
    
def compute_rescaling_values(simulation1, real1):
    """Compute the scaling and alignment values between image1 and mask1"""
    # Flip image1 along the x-axis
    image1_flipped = cv2.flip(simulation1, 0)

    # Find contours in both images
    contours1, _ = cv2.findContours(image1_flipped, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    contours2, _ = cv2.findContours(real1, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    # Compute bounding boxes for both contours
    x1, y1, w1, h1 = cv2.boundingRect(contours1[0])
    x2, y2, w2, h2 = cv2.boundingRect(contours2[0])

    # Compute scaling factors
    scale_x = w2 / w1
    scale_y = h2 / h1
    
    # Uncomment below to visulaise image 1 scaled to image 2 with same axis
    
    '''
    # Rescale image1 using computed scaling factors adn extract only fleece section
    image1_cropped = image1_flipped[y1:y1 + h1, x1:x1 + w1]
    rescaled_image1 = cv2.resize(image1_cropped, (int(image1_cropped.shape[1] * scale_x), int(image1_cropped.shape[0] * scale_y)))
    
    # Crop only masked section of second image
    mask1_cropped = real1[y2:y2 + h2, x2:x2 + w2]

    # Plot the rescaled image1 and mask1 side by side
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    axes[0].imshow(rescaled_image1, cmap='gray')
    axes[0].set_title('Rescaled Image 1')
    axes[0].axis('on')

    axes[1].imshow(mask1_cropped, cmap='gray')
    axes[1].set_title('Original Mask 1')
    axes[1].axis('on')

    plt.tight_layout()
    plt.show()
    
    '''
    return scale_x, scale_y


def rescale_image(image, scale_x, scale_y):
    """Rescales image using the provided scale factors and applies the translation"""
    height, width = image.shape[:2]    
    
    # Rescale image using values from the different between the first images
    scaled_image = cv2.resize(image, (int(width * scale_x), int(height * scale_y)))

    return scaled_image

def compute_percentage_difference(image1, image2):
    """Computes the percentage difference between a specific rectangular area around the center of two binary images"""

    # Define the coordinates where binary comparison will be done
    x_start = 200
    x_end = 600
    y_start = 500
    y_end = 800

    # Crop the specified section from the first image
    cropped_image1 = image1[y_start:y_end, x_start:x_end]
  
    # Remove background from second
    contours, _ = cv2.findContours(image2, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    # Compute bounding boxes for the second image
    x2, y2, w2, h2 = cv2.boundingRect(contours[0])

    # Crop only masked section of second image
    cropped_image2 = image2[y2:y2 + h2, x2:x2 + w2]
    cropped_image2 = cropped_image2[y_start:y_end, x_start:x_end]
    
    # Uncomment below to visualise which section will be compared to one another 
    
    '''
    # Display the cropped regions for both images
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    axes[0].imshow(cropped_image1, cmap='gray')
    axes[0].set_title('Cropped Region from Image 1')
    axes[0].axis('on')

    axes[1].imshow(cropped_image2, cmap='gray')
    axes[1].set_title('Cropped Region from Image 2')
    axes[1].axis('on')

    plt.tight_layout()
    plt.show()
    '''
   
    # Compute absolute difference between the two cropped images
    difference = cv2.absdiff(cropped_image1, cropped_image2)

    # Calculate the number of different pixels
    num_different_pixels = np.sum(difference != 0)

    # Calculate the percentage difference
    total_pixels = np.prod(difference.shape)
    percentage_difference = (num_different_pixels / total_pixels) * 100

    return percentage_difference

def process_images(simulation1, simulation2, real1, real2):
    # Load the images
    simulation1 = load_image(simulation1)
    simulation2 = load_image(simulation2)
    real1 = load_image(real1)
    real2 = load_image(real2)
    
    # Fill holes for the simulation images
    # This is because MPM has some small gaps if not enoguh particles were used in simulation
    simulation1 = fill_holes(simulation1)
    simulation2 = fill_holes(simulation2)
    
    # Compute the rescaling values between image1 and mask1
    rescaling_values = compute_rescaling_values(simulation1, real1)
    if rescaling_values is None:
        print("Unable to compute rescaling values due to missing contours")
        return None

    scale_x, scale_y = rescaling_values
    
    # Apply these values to image2
    scaled_simulation2 = rescale_image(simulation2, scale_x, scale_y)
    
    # Compute the percentage difference between the rescaled image2 and mask2
    percentage_difference = compute_percentage_difference(scaled_simulation2, real2)
    
    return percentage_difference
