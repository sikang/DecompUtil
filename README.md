### MRSL Decomputil Library
Fast convex decomposition on a point cloud. In the basic pipeline, it implements ellipsoid based regional inflation to model free space from a given path.
Detials of the algorithm is proposed in "S. Liu, M. Watterson, K. Mohta, K. Sun, S. Bhattacharya, C.J. Taylor and V. Kumar. Planning Dynamically Feasible Trajectories for Quadrotors using Safe Flight Corridors in 3-D Complex Environments. ICRA 2017".

## Compilation
A) Simple cmake
```sh
$ mkdir build && cd build && cmake .. && make
```

B) Using CATKIN
```sh
$ cd mv decomp_util ~/catkin_ws/src
$ cd ~/catkin_ws & catkin_make -DCMAKE_BUILD_TYPE=Release
```

## Example
![Visualization](./sample/sample1.png)
![Visualization](./sample/sample2.png)

## Doxygen
For more details, please refer to https://sikang.github.io/DecompUtil/index.html

