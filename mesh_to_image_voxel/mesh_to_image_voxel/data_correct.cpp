#include<iostream>
#include<string>
#include<vector>
#include<math.h>
#include<fstream>
#include<time.h>
#include "OpenMeshWarpper.h"
#include "MatVec.h"
#include<Eigen/Dense>
using namespace std;


extern void ergodicfile(const std::string &root, std::vector<std::string> &vecFile, const std::string &fileType);
extern void get_all_face_points_index(Triangle_mesh &mesh, std::vector<OpenMesh::Vec3i> &face_points_idx);


int save_data(Triangle_mesh &mesh, const std::string &save_path)
{

	std::string  save_path_ = save_path + ".stl";
	try
	{
		
		if (!OpenMesh::IO::write_mesh(mesh, save_path_, OpenMesh::IO::Options::Binary))
		{
			std::cerr << "Cannot write mesh to file 'output.obj'" << std::endl;
			return 1;
		}

	}
	catch (std::exception& x)
	{
		std::cerr << x.what() << std::endl;
		return 1;
	}
	return 1;
}


pmp::mat4 AlignTooth(const Eigen::Vector3f & start, const Eigen::Vector3f& end)
{
	float vecdot = start.dot(end);
	Eigen::Vector3f newend = end;
	pmp::mat4 rt = pmp::mat4::identity();
	Eigen::Vector3f rotaxis = start.cross(newend);
	float rotangle = 0.0f;
	if (std::sqrt(rotaxis.dot(rotaxis)) > 0.001f)
	{
		rotaxis = rotaxis.normalized();
		float vecdot = start.dot(newend);
		rotangle = 180.0f*acos(vecdot) / M_PI;
		rt = pmp::rotation_matrix(pmp::vec3(rotaxis[0], rotaxis[1], rotaxis[2]), rotangle);
	}
	return rt;
}


pmp::mat4 get_rotate_matrix(const Eigen::Vector3f & start, const Eigen::Vector3f& end, const Eigen::Vector3f& cpoint)
{

	pmp::mat4 rt1 = pmp::translation_matrix(pmp::vec3(-cpoint[0], -cpoint[1], -cpoint[2]));

	pmp::mat4 rot_mat = AlignTooth(start, end);
	//得到平移矩阵
	pmp::mat4 rt2 = pmp::translation_matrix(pmp::vec3(cpoint[0], cpoint[1], cpoint[2]));

	pmp::mat4 rtmatrix = rt2 * rot_mat* rt1;

	return rtmatrix;
}

Triangle_mesh mesh_rotate_transform(Triangle_mesh &mesh, const pmp::mat4 &rot)
{
	Triangle_mesh::Point cen(0, 0, 0);
#pragma omp parallel for
	for (auto vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
	{
		Triangle_mesh::Point point = mesh.point(*vi);
		pmp::vec4 pv4_point = pmp::vec4(point[0], point[1], point[2], 1);
		pv4_point = rot * pv4_point;
		mesh.set_point(*vi, Triangle_mesh::Point(pv4_point[0], pv4_point[1], pv4_point[2]));
		cen = cen + point;
	}
	std::cout << cen / mesh.n_vertices() << std::endl;

	std::string save_path = "correct_data";
	save_data(mesh, save_path);


	return mesh;
}





void cal_mesh_averge_normal(Triangle_mesh &mesh, Eigen::Vector3f &aver_normal)
{
	//得到所有face的点的顶点id号
	std::vector<OpenMesh::Vec3i> face_points_idx;
	get_all_face_points_index(mesh, face_points_idx);

	int vertic_nums = mesh.n_vertices();
	int face_nums = mesh.n_faces();
	auto points = mesh.points();
	std::vector<OpenMesh::Vec3f> all_points;
	for (int i = 0; i < vertic_nums; i++)
	{
		all_points.push_back(points[i]);
	}

	OpenMesh::Vec3f normal_(0, 0, 0);
	OpenMesh::Vec3i fv_idx(0, 0, 0);
	OpenMesh::Vec3f v0(0, 0, 0), v1(0, 0, 0), v2(0, 0, 0);
	for (int fi = 0; fi < face_points_idx.size(); fi++)
	{
		fv_idx = face_points_idx[fi];
		v0 = all_points[fv_idx[0]];
		v1 = all_points[fv_idx[1]];
		v2 = all_points[fv_idx[2]];

		//求解平面方程，使用三点，返回法向量//cal normal
		OpenMesh::Vec3f vec1 = OpenMesh::Vec3f(v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]);
		OpenMesh::Vec3f vec2 = OpenMesh::Vec3f(v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]);
		normal_ = normal_ + vec1.cross(vec2).normalized();
	}
	normal_ = normal_.normalized();
	aver_normal = Eigen::Vector3f(normal_[0], normal_[1], normal_[2]);
}


void new_data_coorrect(Triangle_mesh &mesh, pmp::mat4 &rot_mat)
{
	int point_nums = mesh.n_vertices();
	Triangle_mesh::Point cpoint(0, 0, 0);
	for (auto v_irer = mesh.vertices_begin(); v_irer != mesh.vertices_end(); ++v_irer)
	{
		cpoint = cpoint + mesh.point(*v_irer);
	}
	//center point
	Eigen::Vector3f cenpoint(cpoint[0] / point_nums, cpoint[1] / point_nums, cpoint[2] / point_nums);

	Eigen::Vector3f aver_normal(0, 0, 1);
	cal_mesh_averge_normal(mesh, aver_normal);


	Eigen::Vector3f end(0, 0, 1);
	rot_mat = get_rotate_matrix(aver_normal, end, cenpoint);

	Triangle_mesh correct_mesh = mesh_rotate_transform(mesh, rot_mat);



}



int main(int argc, char *argv[])
{


	std::string fileType = ".stl";
	std::vector<std::string> vecFile;
	const std::string data_root = "./";
	std::string save_root = "./";
	ergodicfile(data_root, vecFile, fileType);

	for (int i = 0; i < 1; i++)
	{
		vecFile[i] = "./data.stl";

		Triangle_mesh mesh;

		if (!OpenMesh::IO::read_mesh(mesh, vecFile[i]))
		{
			return 0;
		}
		//数据矫正
		pmp::mat4 op_rot = pmp::mat4::identity();
		new_data_coorrect(mesh, op_rot);

		std::cout << vecFile[i] << std::endl;
	}

	return 0;
}
