#include "functions.h"
#include <iostream>
#include <fstream>
#include <vector>


// Define TM as a 3x3 matrix of floats
struct TM {
    float data[3][3];
};

// Define dim based on the dimensionality of your data
const int dim = 3;

void readF(const std::string& filename, std::vector<TV>& Xp_data) {
    std::ifstream infile(filename, std::ios::binary | std::ios::in);
    if (!infile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    int nump;
    infile.read((char*)&nump, sizeof(int));

    // Resize the Xp_data vector based on nump
    Xp_data.resize(nump);
    
    std::vector<TM> def_data(nump); // Temporary storage for deformation data, can be removed if not needed

    for (int i = 0; i < nump; i++) {
        TV Xp;
        for (int q = 0; q < dim; q++) {
            float temp;
            infile.read((char*)&temp, sizeof(float));
            if (q == 0) Xp.x = temp;
            else if (q == 1) Xp.y = temp;
            else Xp.z = temp;
        }
        
        def_data[i].data[0][0] = def_data[i].data[0][1] = def_data[i].data[0][2] =
            def_data[i].data[1][0] = def_data[i].data[1][1] = def_data[i].data[1][2] =
            def_data[i].data[2][0] = def_data[i].data[2][1] = def_data[i].data[2][2] = 0.0f;

        for (int p = 0; p < dim; p++) {
            for (int q = 0; q < dim; q++) {
                float temp;
                infile.read((char*)&temp, sizeof(float));
                def_data[i].data[q][p] = temp;
            }
        }

        Xp_data[i] = Xp;
    }

    infile.close();

}
