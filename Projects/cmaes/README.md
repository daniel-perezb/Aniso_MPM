## TODO
- [ ] Add code for obtaining fibre directons using structure tensors
- [ ] Add code to obtain ply file for simulation using first frame of video and mask
- [ ] How is the output from TAPIR going to be rescaled now
- [ ] Address all the problems in image comparison


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
To run the training of the parameters run from the 'cmaes' folder

    python3 cmaes_animation.py 

## Image Comparison
Image comparison is done by creating an image and rescaling from the output files named 'data_*.dat'
TODO: 
- [ ] Change cmaes_animation.py to perform image comparison automaticaly. 
- [ ] Scale values needed for MSE must be obtained automatically

## Visualization
An image of the last pointcloud is displayed using main_image.cpp in calculater_mse. This image is then stored as 'binary_image.png'. 
THe rest of the files are '.bgeo' files compatible with Houdini's VDB software. 


