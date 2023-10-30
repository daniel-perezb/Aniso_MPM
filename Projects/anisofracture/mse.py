import numpy as np
import os

from read_data import readF
from scale_data import readData, readScaleFactors, rescaleAndTranslateData


def distance2D(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def findClosestPoint(interestPoint, XP_data):
    minDistance = float('inf')
    closestPoint = None
    index = -1

    for i, point in enumerate(XP_data):
        distance = distance2D(interestPoint[0], interestPoint[1], point[0], point[2])
        if distance < minDistance:
            minDistance = distance
            closestPoint = point
            index = i

    return closestPoint, index


def calculate_mse(dataFilePath, logicFilePath, scaleFilePath, directory):
    # Reading data from files
    data = readData(dataFilePath)
    logic_data = readData(logicFilePath)
    scaleFactors = readScaleFactors(scaleFilePath)
    rescaledData = rescaleAndTranslateData(data, scaleFactors)

    Xp_data = readF(f"{directory}/data_0.dat")

    # Find the closest point for each point in rescaledData
    indexOfClosestPoints = []
    all_data = []

    for row in rescaledData:
        if len(row) >= 2:
            interestPoint = [row[0], row[1]]
            closestPoint, indexOfClosestPoint = findClosestPoint(interestPoint, Xp_data)
            indexOfClosestPoints.append(indexOfClosestPoint)
            all_data.append(closestPoint)

    # List all files in the directory
    all_files = os.listdir(directory)

    # Filter out files that match the "data_{}.dat" pattern and extract numbers
    file_numbers = [
        int(f.split("_")[1].split(".")[0])
        for f in all_files
        if f.startswith("data_") and f.endswith(".dat")
    ]

    # Get the highest number
    highest_number = max(file_numbers)
    filename = f"{directory}/data_{highest_number}.dat"
    currentXp_data = readF(filename)
    video_data = []
    for index in indexOfClosestPoints:
        if 0 <= index < len(currentXp_data):
            video_data.append(currentXp_data[index])

    # Extract the data from the last frame
    last_frame_data = []
    for row in rescaledData:
        if len(row) > highest_number:  # Ensure the row has enough elements
            interestPoint = [row[highest_number], row[highest_number + 1]]  # Extract pair
            last_frame_data.append(interestPoint)

    # Convert to numpy arrays
    last_frame_data = np.array(last_frame_data)
    video_data = np.array(video_data)

    # Extract first and last column video data
    video_data = video_data[:, [0, 2]]

    # Calculate MSE
    mse = np.square(np.subtract(last_frame_data, video_data)).mean()

    return mse
