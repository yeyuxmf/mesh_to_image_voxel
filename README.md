# mesh_to_image_voxel
Mesh to depth image, Mesh to rendering image, Mesh to voxel，3D voxel reconstruction，Eigen, OpenMesh, C++, Python

# First step
1. Firstly, it is necessary to perform directional correction on the original mesh data. Ensure that the sum of the normals of all triangular surfaces faces upwards.
2. Visualize as follows:  data_correct.cpp   
![correct](https://github.com/huang229/mesh_to_image_voxel/assets/29627190/2d9ea297-a37e-4674-8a4f-c7fe5f8db84a)

# Mesh to depth image, Mesh to rendering image 
1. References：https://zhuanlan.zhihu.com/p/363245957
2.  mesh_to_image_voxel.cpp  mesh_to_image(vecFile[i], save_root)
![image](https://github.com/huang229/mesh_to_image_voxel/assets/29627190/8686a642-9c30-4196-ba1b-17730ab1fdb2)

# Mesh to voxel  
1. 3Dmesh voxelization only results in one layer of shell. If you want to obtain a solid voxelization result, you need to modify the code slightly.       
2. mesh_to_image_voxel.cpp mesh_transform_convex_hull(vecFile[i], save_root)    
![voxel](https://github.com/huang229/mesh_to_image_voxel/assets/29627190/25b5760e-600e-43f8-a425-11c4c6e16c87)






