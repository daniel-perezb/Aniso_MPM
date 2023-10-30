import struct

def readF(filename):
    Xp_data = []
    with open(filename, 'rb') as infile:
        nump = struct.unpack('i', infile.read(4))[0]  # Read integer

        for _ in range(nump):
            Xp = []
            for _ in range(3):  # 3 dimensions: x, y, z
                temp_float = struct.unpack('f', infile.read(4))[0]
                Xp.append(temp_float)
            
            # If the deformation data (TM) is important, add reading code for it here.
            # For now, I'm skipping it based on your initial conversion.

            Xp_data.append(Xp)

    return Xp_data