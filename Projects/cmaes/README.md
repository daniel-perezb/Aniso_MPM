## TODO
- [ ] Check if code cimpiles in a different computer
- [ ] Upload parameters used for paper


## Code description 

This is the opensource code for the following works:

(1) Non-rigid image registration for elastoplastic materials using physics-based simulation, Daniel Perez, Raphael Falque, Alen Alempijevic 

The paper above utilises the material simulation implemented below for obtaining material's properties of fleece:

(2) AnisoMPM: Animating Anisotropic Damage Mechanics, Joshuah Wolper, Yunuo Chen, Minchen Li, Yu Fang, Ziyin Qu, Jiecong Lu, Meggie Cheng, Chenfanfu Jiang (SIGGRAPH 2020)
Project Page: https://joshuahwolper.com/anisompm

## Data
Data for the animation can be downloaded in the following link: 
https://drive.google.com/file/d/1U5400w3e5UFbbBOVgrsYPe5iXhBnDi4C/view?usp=drive_link

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
To perform pixel si ilarity between the output and the target Houdini software must be installed. 
Hyperparameters in lines 116-119 in 'cmaes-animation' must be manually adjusted.

If wanting to run the animation without the pixel similarity keep variable equal to false in line 13 of 'cmaes_animation', otherwise:
- [ ] Explain how to adjust them 

## Visualization
All files are output as Houdini's VDB '.bgeo' files. 


