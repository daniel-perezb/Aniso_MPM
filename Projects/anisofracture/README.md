## Introduction

Welcome to AnisoMPM!

This code distribution contains one 2D demo (diskTear) and 13 3D demos (listed below):
1. diskShoot
2. cheese
3. parameterEffects
4. dongpoPork
5. hangingTorus
6. implicitVsExplicit
7. brokenHeart
8. boneTwist
9. tubeCompress
10. diskTear
11. tubePull
12. orange
13. fish
14. fleece

## Instructions for Use

0. Download (and unzip in root directory of ziran2020) all of the data directory found here: https://www.seas.upenn.edu/~cffjiang/research/ziran2020/Data.zip

1. To switch between 2D/3D demos, go into anisofracture.cpp and change line 14 to define the problem dimension.

2. Open anisofractureBatch2D.py or anisofractureBatch3D.py and set the controls within them to dictate which demos will be run.

3. Make sure all files are downloaded and correctly organized in either Data/TetMesh (for .mesh) or Data/LevelSets (for .vdb)

4. Compile in build directory as "make anisofracture -j8"

5. To run your selected demos, run "python anisofractureBatch3D.py" (for 3D).


**Each demo can also be changed directly in the examples folder, but the Python approach requires minimal
compile time: simply compile once and all demos should be runnable through the Python script! 


