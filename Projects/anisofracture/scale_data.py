def readData(filePath):
    data = []
    with open(filePath, "r") as file:
        for line in file:
            currentLineData = [float(x) for x in line.split()]
            data.append(currentLineData)
    return data

def readScaleFactors(filePath):
    with open(filePath, "r") as file:
        factors = file.readline().split(',')
    return factors

def rescaleAndTranslateData(data, scaleFactors):
    if len(scaleFactors) != 6:
        print("Invalid number of scale factors.")
        return data

    scaleX, _, scaleZ, translateX, _, translateZ = map(float, scaleFactors)

    result = []
    for row in data:
        newRow = []
        for i in range(0, len(row), 2):
            newRow.append((row[i] - translateX) / scaleX)
            newRow.append((row[i+1] - translateZ) / scaleZ)
        result.append(newRow)
    return result
