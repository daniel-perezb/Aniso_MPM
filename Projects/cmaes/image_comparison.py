import cv2
import numpy as np

def load_image(image_path):
    """Load an image in grayscale."""
    return cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)

def compute_rescaling_values(image1, mask1):
    """Compute the scaling and alignment values between image1 and mask1."""
    # Find contours in both images
    contours1, _ = cv2.findContours(image1, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    contours2, _ = cv2.findContours(mask1, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    # Compute bounding boxes for both contours
    x1, y1, w1, h1 = cv2.boundingRect(contours1[0])
    x2, y2, w2, h2 = cv2.boundingRect(contours2[0])
    
    # Compute scaling factors
    scale_x = w2 / w1
    scale_y = h2 / h1
    
    return scale_x, scale_y, x1, y1, x2, y2

def rescale_image(image, scale_x, scale_y, shift_x, shift_y):
    """Rescales image using the provided scale factors and applies the translation."""
    height, width = image.shape[:2]
    
    # Apply scaling
    scaled_image = cv2.resize(image, (int(width * scale_x), int(height * scale_y)))
    
    # Apply translation
    translation_matrix = np.float32([[1, 0, shift_x], [0, 1, shift_y]])
    translated_image = cv2.warpAffine(scaled_image, translation_matrix, (width, height))
    
    return translated_image

def compute_percentage_difference(image1, image2):
    """Computes the percentage difference between two binary images."""
    # Compute absolute difference between the two images
    difference = cv2.absdiff(image1, image2)
    
    # Calculate the number of different pixels
    num_different_pixels = np.sum(difference != 0)
    
    # Calculate the percentage difference
    total_pixels = np.prod(difference.shape)
    percentage_difference = (num_different_pixels / total_pixels) * 100
    
    return percentage_difference

def process_images(image1, image2, mask1, mask2):

    # Load the images
    image1 = load_image(image1)
    image2 = load_image(image2)
    mask1 = load_image(mask1)
    mask2 = load_image(mask2)
    
    #Compute the rescaling values between image1 and mask1
    scale_x, scale_y, x1, y1, x2, y2 = compute_rescaling_values(image1, mask1)
    
    #Apply these values to image2
    rescaled_image2 = rescale_image(image2, scale_x, scale_y, x2 - x1, y2 - y1)
    
    #Compute the percentage difference between the rescaled image2 and mask2
    percentage_difference = compute_percentage_difference(rescaled_image2, mask2)
    
    return percentage_difference
