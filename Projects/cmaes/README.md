## Code description 

This is the opensource code for the following works:

(1) 

The paper above utilises the material simulation implemented below for obtaining material's properties of fleece:

(2) AnisoMPM: Animating Anisotropic Damage Mechanics, Joshuah Wolper, Yunuo Chen, Minchen Li, Yu Fang, Ziyin Qu, Jiecong Lu, Meggie Cheng, Chenfanfu Jiang (SIGGRAPH 2020)
Project Page: https://joshuahwolper.com/anisompm

## Data
Data for the animation can be downloaded in the following link: 
https://drive.google.com/file/d/1Cfn29vXSrPEiAoZ3FlQq0cbJG4EshmSd/view?usp=drive_link

The data must be unzip in the following directory: 
/Fibre_Directions/Aniso_MPM/Data/TetMesh/

## Path change
Ensure that line 56 is set correctly in file '../Projects/anisofracture/examples/fleece.h
- [ ] Run training to check if code is working correctly

## Single Run
To run the code only once:
    
    cd ../anisofracture
    ./anisofracture -test 27
   

Material properties of the fleece are located in 'parameters.txt' 


## Dependencies Installation

    # Update the package list
    sudo apt-get update

    # Install Python and development tools
    sudo apt-get install -y python3-pip python3-dev build-essential libopencv-dev 

    # Install Python packages   
    pip install numpy cmaes opencv-python

   
## Building MSE

    cd calculate_mse
    mkdir build && cd build
    cmake .. && make -j 4

## Running training of cmaes
To run the training of the parameters run from the 'cmaes' folder.

    python3 cmaes_animation.py 


## Creating mesh from image
To create a mesh from an image the file 'mesh_creation.py' is used. This files takes in a binary image of the object and the camera extrinsics to create a mesh with real-life dimensions. THe output file is then transformed into a '*.mesh' file by using TetWild library

https://github.com/Yixin-Hu/TetWild


## Fibre Direction
Fibre direction for fleece images is perfomed using structure tensors. Direction of the fibres can be obtained using 'direction_fibres.py'

## Image Comparison
Image comparison is done by creating an image and rescaling from the output files named 'data_*.dat. A new image is created from running:
    ./calculate_mse/build/create_image 

The images are then stored as output_binary_image1.png and output_binary_image1.png in folder '../../Data/TetMesh/fleece_files/'

## Visualization
An image of the last pointcloud is displayed using main_image.cpp in calculate_mse. This image is then stored as previously discussed. 
The rest of the files are '.bgeo' files compatible with Houdini's VDB software. 


