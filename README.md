# mesh_to_image_voxel
    Mesh to depth image, Mesh to rendering image, Mesh to voxel，3D voxel reconstruction，Eigen, OpenMesh, C++, Python
1. In addition to the extremely slow speed of 3D voxel reconstruction, other speeds are extremely fast, and you can further optimize according to your own needs.

# Visualization Tool Installation Package
链接：https://pan.baidu.com/s/1SMPms9_GoujxcmbTCV0K3Q  （Geomagic Wrap 2013 (64 bit)）
提取码：fiid  
This is a 3D visualization tool that supports visualization of stl/obj files and point cloud txt file formats.
# Environment
1. vs2017
2. eigen-3.4.0
3. OpenMesh 9.1

# First step
1. Firstly, it is necessary to perform directional correction on the original mesh data. Ensure that the sum of the normals of all triangular surfaces faces upwards.
2. Visualize as follows:  data_correct.cpp   
![correct](https://github.com/huang229/mesh_to_image_voxel/assets/29627190/2d9ea297-a37e-4674-8a4f-c7fe5f8db84a)

# Mesh to depth image, Mesh to rendering image 
1. References：https://zhuanlan.zhihu.com/p/363245957
2.  mesh_to_image_voxel.cpp  mesh_to_image(vecFile[i], save_root)
3.  Visualization code: Under the visualization_python directory
![image](https://github.com/huang229/mesh_to_image_voxel/assets/29627190/8686a642-9c30-4196-ba1b-17730ab1fdb2)

# Mesh to voxel  
1. 3Dmesh voxelization only results in one layer of shell. If you want to obtain a solid voxelization result, you need to modify the code slightly.       
2. mesh_to_image_voxel.cpp mesh_transform_convex_hull(vecFile[i], save_root)    
![voxel](https://github.com/huang229/mesh_to_image_voxel/assets/29627190/25b5760e-600e-43f8-a425-11c4c6e16c87)


# 3D voxel reconstruction
1. This implementation mainly refers to the existing code on the internet, but the source has been forgotten. If there are any licensing rights, please follow the original author's permission.
2. Below is a 3D reconstruction of head scan data from CBCT.
![reconstruction](https://github.com/huang229/mesh_to_image_voxel/assets/29627190/ef84a62c-4baf-48c4-953c-f6726cf4c4ce)


# License and Citation
1.Without permission, the design concept of this model shall not be used for commercial purposes, profit seeking, etc.

2.If you refer to the design concept of this model for theoretical research, please also add a reference.





