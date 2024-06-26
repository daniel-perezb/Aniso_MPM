import cv2
import numpy as np


def preprocess_image(image_path, scale=1.0, translation=(0, 0)):
    # Load the image
    image = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)
    if image is None:
        raise ValueError(f"Image not found at path: {image_path}")

    # Apply scaling and translation
    if scale != 1.0:
        width = int(image.shape[1] * scale)
        height = int(image.shape[0] * scale)
        image = cv2.resize(image, (width, height), interpolation=cv2.INTER_AREA)
    if any(translation):
        M = np.float32([[1, 0, translation[0]], [0, 1, translation[1]]])
        image = cv2.warpAffine(image, M, (image.shape[1], image.shape[0]))

    # Convert to grayscale and apply threshold
    gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    _, binary_image = cv2.threshold(gray_image, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
    return binary_image


def pad_to_match(img1, img2):
    # Pad the smaller image to match the size of the larger one
    h1, w1 = img1.shape[:2]
    h2, w2 = img2.shape[:2]

    # Determine padding amounts
    vertical_padding = abs(h1 - h2)
    horizontal_padding = abs(w1 - w2)

    # Pad the smaller image
    if h1 > h2:
        img2 = cv2.copyMakeBorder(img2, vertical_padding // 2, vertical_padding - vertical_padding // 2,
                                  0, 0, cv2.BORDER_CONSTANT, value=[0, 0, 0])
    elif h2 > h1:
        img1 = cv2.copyMakeBorder(img1, vertical_padding // 2, vertical_padding - vertical_padding // 2,
                                  0, 0, cv2.BORDER_CONSTANT, value=[0, 0, 0])

    if w1 > w2:
        img2 = cv2.copyMakeBorder(img2, 0, 0,
                                  horizontal_padding // 2, horizontal_padding - horizontal_padding // 2,
                                  cv2.BORDER_CONSTANT, value=[0, 0, 0])
    elif w2 > w1:
        img1 = cv2.copyMakeBorder(img1, 0, 0,
                                  horizontal_padding // 2, horizontal_padding - horizontal_padding // 2,
                                  cv2.BORDER_CONSTANT, value=[0, 0, 0])

    return img1, img2
def compute_percentage_difference(target_image_path, scale, translation, rect_coords, screenshot_coords):
    # Load the image
    comparison_image_path = '/home/daniel/Fleece_animations/Extract_STL/screenshot.png'

    # Preprocess images
    binary_target_image = preprocess_image(target_image_path)
    binary_comparison_image = preprocess_image(comparison_image_path, scale=scale, translation=translation)

    # Pad images to match sizes
    binary_target_image, binary_comparison_image = pad_to_match(binary_target_image, binary_comparison_image)

    # Extract the ROI based on rect_coords
    x, y, w, h = rect_coords
    roi_target = binary_target_image[y:y + h, x:x + w]
    roi_comparison = binary_comparison_image[y:y + h, x:x + w]

    # Compute the pixel difference within the rectangle
    difference = cv2.absdiff(roi_target, roi_comparison)
    num_different_pixels = np.sum(difference != 0)

    # Calculate the percentage difference
    total_pixels = np.prod(difference.shape)
    percentage_difference = (num_different_pixels / total_pixels) * 100

    return percentage_difference



