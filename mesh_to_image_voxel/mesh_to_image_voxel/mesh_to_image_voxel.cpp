#include<iostream>
#include<string>
#include<vector>
#include<math.h>
#include<fstream>
#include<time.h>
#include "OpenMeshWarpper.h"
//#include <Utils/MatVec.h>
using namespace std;


extern void ergodicfile(const std::string &root, std::vector<std::string> &vecFile, const std::string &fileType);

void get_all_face_points_index(Triangle_mesh &mesh, std::vector<OpenMesh::Vec3i> &face_points_idx)
{
	int vertic_nums = mesh.n_vertices();
	int face_nums = mesh.n_faces();

#pragma omp parallel for
	for (auto f : mesh.faces()) {

		auto fv = mesh.cfv_begin(f);

		int index0 = fv->idx(); ++fv;
		int index1 = fv->idx(); ++fv;
		int index2 = fv->idx();
		face_points_idx.push_back(OpenMesh::Vec3i(index0, index1, index2));
	}


}

void transform_points(std::vector<OpenMesh::Vec3f> &all_points, const Triangle_mesh::Point &min_point)
{
#pragma omp parallel for
	for (int iv = 0; iv < all_points.size(); iv++)
	{
		//平移0.5个单位，离边界一定距离
		all_points[iv] = all_points[iv] - min_point + Triangle_mesh::Point(0.5, 0.5, 0.5);
	}

}

void transform_center(std::vector<OpenMesh::Vec3f> &all_points, const Triangle_mesh::Point &max_point,
	const Triangle_mesh::Point &min_point, Triangle_mesh::Point &trans_point)
{
	Triangle_mesh::Point length = max_point - min_point;
	float max_length = length.max();
#pragma omp parallel for
	for (int iv = 0; iv < all_points.size(); iv++)
	{
		OpenMesh::Vec3f vp = all_points[iv];
		vp[0] = vp[0] + (max_length - length[0]) / 2.0;
		vp[1] = vp[1] + (max_length - length[1]) / 2.0;
		vp[2] = vp[2] + (max_length - length[2]) / 2.0;
		all_points[iv] = vp;

	}
	trans_point[0] = (max_length - length[0]) / 2.0;
	trans_point[1] = (max_length - length[1]) / 2.0;
	trans_point[2] = (max_length - length[2]) / 2.0;

}



//如果三角形ABP、BCP、CAP都是按照逆时针的方向排列的，那么P必然在ABC的外部。
//(x1,y1,0) X (x2,y2,0) = (0, 0, x1y2 – y1x2)  添加一个维度，使用叉乘判断是否在三角形内部
bool IsInTriangle(float x, float y, float ax, float ay, float bx, float by, float cx, float cy)
{
	return  (bx - ax) * (y - ay) > (by - ay) * (x - ax) &&
		(cx - bx) * (y - by) > (cy - by) * (x - bx) &&
		(ax - cx) * (y - cy) > (ay - cy) * (x - cx) ? false : true;

}

bool IsInTriangle_(float x, float y, float x1, float y1, float x2, float y2, float x3, float y3)
{

	float y2_y3 = y2 - y3;
	float x3_x2 = x3 - x2;
	float x1_x3 = x1 - x3;
	float y1_y3 = y1 - y3;
	float y_y3 = y - y3;
	float x_x3 = x - x3;
	float y3_y1 = y3 - y1;

	float fmu = y2_y3 * x1_x3 + x3_x2 * y1_y3;

	float nmn1 = (y2_y3 * x_x3 + x3_x2 * y_y3) / fmu;
	float nmn2 = (y3_y1 * x_x3 + x1_x3 * y_y3) / fmu;

	float nmn3 = 1 - nmn1 - nmn2;

	if (nmn1 >= 0 && nmn1 <= 1 && nmn2 >= 0 && nmn2 <= 1 && nmn3 >= 0 && nmn3 <= 1)
	{
		return true;
	}
	else
	{
		return false;
	}

}


void mesh_project_xyplane(const std::vector<OpenMesh::Vec3i> &face_points_idx, const std::vector<OpenMesh::Vec3f> &all_points,
	float* map_img, float* render_img, const float &resolution, const int &width, const int &height)
{
	OpenMesh::Vec3i fv_idx(0, 0, 0);
	OpenMesh::Vec3f v0(0, 0, 0), v1(0, 0, 0), v2(0, 0, 0);
	OpenMesh::Vec2f pv0(0, 0), pv1(0, 0), pv2(0, 0);
	OpenMesh::Vec3f dir = OpenMesh::Vec3f(0, 0, 1);
	for (int fi = 0; fi < face_points_idx.size(); fi++)
	{
		fv_idx = face_points_idx[fi];
		v0 = all_points[fv_idx[0]];
		v1 = all_points[fv_idx[1]];
		v2 = all_points[fv_idx[2]];

		//求解平面方程，使用三点，返回法向量
		OpenMesh::Vec3f vec1 = OpenMesh::Vec3f(v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]);
		OpenMesh::Vec3f vec2 = OpenMesh::Vec3f(v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]);
		OpenMesh::Vec3f normal_ = vec1.cross(vec2).normalized();

		//投影到xy平面，在标准世界坐标系中，那么就是忽略Z轴值,得到平面三角形
		//并反算为平面坐标范围内的值
		pv0 = OpenMesh::Vec2f(v0[0], v0[1]).operator*(resolution);
		pv1 = OpenMesh::Vec2f(v1[0], v1[1]).operator*(resolution);
		pv2 = OpenMesh::Vec2f(v2[0], v2[1]).operator*(resolution);
		v0 = v0.operator*(resolution);


		OpenMesh::Vec2f min_v = pv0;
		OpenMesh::Vec2f max_v = pv0;
		min_v = min_v.minimize(pv1);
		min_v = min_v.minimize(pv2);
		max_v = max_v.maximize(pv1);
		max_v = max_v.maximize(pv2);

		//得到三角形在xy平面张成的长方形范围值

		int w_st = min_v[0];
		int w_ed = max_v[0];
		int h_st = min_v[1];
		int h_ed = max_v[1];


		for (int y = h_st; y <= h_ed; y++)
		{
			for (int x = w_st; x <= w_ed; x++)
			{

				//如果三角形ABP、BCP、CAP都是按照逆时针的方向排列的，那么P必然在ABC的外部
				bool is_flag = IsInTriangle_(x, y, pv0[0], pv0[1], pv1[0], pv1[1], pv2[0], pv2[1]);
				if (true == is_flag)
				{
					//v0
					float x_x0 = x - v0[0];
					float y_y0 = y - v0[1];
					float value = normal_[0] * x_x0 + normal_[1] * y_y0 + normal_[2] * v0[2];
					float z_value = value / normal_[2]; //这里没考虑面片垂直xy平面

					int index_xy = y * width + x;
					if (map_img[index_xy] < z_value)
					{
						map_img[index_xy] = z_value;
						float colorv = dir.dot(normal_);
						render_img[index_xy] = colorv > 0 ? colorv : 0;
					}

				}


			}

		}

	}

}




void mesh_voxel_shell_(const std::vector<OpenMesh::Vec3i> &face_points_idx, const std::vector<OpenMesh::Vec3f> &all_points,
	float* map_img, const float &resolution, const int &width, const int &height, const int &depth)
{
	OpenMesh::Vec3i fv_idx(0, 0, 0);
	OpenMesh::Vec3f v0(0, 0, 0), v1(0, 0, 0), v2(0, 0, 0), cp(0, 0, 0);
	OpenMesh::Vec3f pv0(0, 0, 0), pv1(0, 0, 0), pv2(0, 0, 0);
	for (int fi = 0; fi < face_points_idx.size(); fi++)
	{
		fv_idx = face_points_idx[fi];
		v0 = all_points[fv_idx[0]];
		v1 = all_points[fv_idx[1]];
		v2 = all_points[fv_idx[2]];

		//求解平面方程，使用三点，返回法向量
		OpenMesh::Vec3f vec1 = OpenMesh::Vec3f(v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]);
		OpenMesh::Vec3f vec2 = OpenMesh::Vec3f(v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]);
		OpenMesh::Vec3f normal_ = vec1.cross(vec2).normalized();

		//投影到xy平面，在标准世界坐标系中，那么就是忽略Z轴值,得到平面三角形
		//并反算为平面坐标范围内的值
		pv0 = OpenMesh::Vec3f(v0[0], v0[1], v0[2]).operator*(resolution);
		pv1 = OpenMesh::Vec3f(v1[0], v1[1], v1[2]).operator*(resolution);
		pv2 = OpenMesh::Vec3f(v2[0], v2[1], v2[2]).operator*(resolution);
		v0 = v0.operator*(resolution);

		cp = (pv0 + pv1 + pv2) / 3;

		OpenMesh::Vec3f min_v = pv0;
		OpenMesh::Vec3f max_v = pv0;
		min_v = min_v.minimize(pv1);
		min_v = min_v.minimize(pv2);
		max_v = max_v.maximize(pv1);
		max_v = max_v.maximize(pv2);

		//得到三角形在xy平面张成的长方形范围值

		int w_st = floor(min_v[0]) > 0 ? floor(min_v[0]) : 0;
		int w_ed = ceil(max_v[0]) < width ? ceil(max_v[0]) : width;
		int h_st = floor(min_v[1]) > 0 ? floor(min_v[1]) : 0;
		int h_ed = ceil(max_v[1]) < height ? ceil(max_v[1]) : height;
		int d_st = floor(min_v[2]) > 0 ? floor(min_v[2]) : 0;
		int d_ed = ceil(max_v[2]) < depth ? ceil(max_v[2]) : depth;

		//yx
		for (int z = d_st; z < d_ed; z++)
		{
			for (int y = h_st; y < h_ed; y++)
			{
				for (int x = w_st; x < w_ed; x++)
				{

					int z0 = z;
					int z1 = z + 1;
					int y0 = y;
					int y1 = y + 1;
					int x0 = x;
					int x1 = x + 1;

					OpenMesh::Vec3f vec0 = OpenMesh::Vec3f(x0 - cp[0], y0 - cp[1], z0 - cp[2]);
					OpenMesh::Vec3f vec1 = OpenMesh::Vec3f(x0 - cp[0], y0 - cp[1], z1 - cp[2]);
					OpenMesh::Vec3f vec2 = OpenMesh::Vec3f(x0 - cp[0], y1 - cp[1], z0 - cp[2]);
					OpenMesh::Vec3f vec3 = OpenMesh::Vec3f(x0 - cp[0], y1 - cp[1], z1 - cp[2]);
					OpenMesh::Vec3f vec4 = OpenMesh::Vec3f(x1 - cp[0], y0 - cp[1], z0 - cp[2]);
					OpenMesh::Vec3f vec5 = OpenMesh::Vec3f(x1 - cp[0], y0 - cp[1], z1 - cp[2]);
					OpenMesh::Vec3f vec6 = OpenMesh::Vec3f(x1 - cp[0], y1 - cp[1], z0 - cp[2]);
					OpenMesh::Vec3f vec7 = OpenMesh::Vec3f(x1 - cp[0], y1 - cp[1], z1 - cp[2]);


					float len0 = vec0.dot(normal_);
					float len1 = vec1.dot(normal_);
					float len2 = vec2.dot(normal_);
					float len3 = vec3.dot(normal_);
					float len4 = vec4.dot(normal_);
					float len5 = vec5.dot(normal_);
					float len6 = vec6.dot(normal_);
					float len7 = vec7.dot(normal_);

					bool sgn0 = (len0 >= 0);
					bool sgn1 = (len1 >= 0);
					bool sgn2 = (len2 >= 0);
					bool sgn3 = (len3 >= 0);
					bool sgn4 = (len4 >= 0);
					bool sgn5 = (len5 >= 0);
					bool sgn6 = (len6 >= 0);
					bool sgn7 = (len7 >= 0);

					bool adj0 = (abs(len0) > 10e-5);
					bool adj1 = (abs(len1) > 10e-5);
					bool adj2 = (abs(len2) > 10e-5);
					bool adj3 = (abs(len3) > 10e-5);
					bool adj4 = (abs(len4) > 10e-5);
					bool adj5 = (abs(len5) > 10e-5);
					bool adj6 = (abs(len6) > 10e-5);
					bool adj7 = (abs(len7) > 10e-5);

					int sumva = adj0 + adj1 + adj2 + adj3 + adj4 + adj5 + adj6 + adj7;
					int sumv = sgn0 + sgn1 + sgn2 + sgn3 + sgn4 + sgn5 + sgn6 + sgn7;

					if (sumv < 8 && sumv >10e-3 && (sumva == 8))
					{
						int index_xy = z * height * width + y * width + x;
						map_img[index_xy] = 1;
					}
				}

			}
		}



	}
}


int mesh_to_image(std::string &file_path, std::string &save_root)
{
	Triangle_mesh mesh;
	bool flag_io = OpenMesh::IO::read_mesh(mesh, file_path);
	if (false == flag_io)
	{
		std::cout << "load data eror" << std::endl;
		return 0;
	}

	double beginTime = clock();
	int vertic_nums = mesh.n_vertices();
	int face_nums = mesh.n_faces();
	auto points = mesh.points();

	Triangle_mesh::Point min_point(1000, 1000, 1000), max_point(-1000, -1000, -1000);

	std::vector<OpenMesh::Vec3f> all_points;
	for (int i = 0; i < vertic_nums; i++)
	{
		min_point = min_point.minimize(points[i]);
		max_point = max_point.maximize(points[i]);
		all_points.push_back(points[i]);
	}

	//得到所有face的点的顶点id号
	std::vector<OpenMesh::Vec3i> face_points_idx;
	get_all_face_points_index(mesh, face_points_idx);

	//减去最小值，使最小点为0
	transform_points(all_points, min_point);


	//投影到固定尺寸二维平面图像，长宽缩放比一致，这里考虑z轴是垂直空间平面xy轴
	//投影还得考虑中心坐标对齐
	int width = 1024;
	int height = 1024;
	Triangle_mesh::Point length = max_point - min_point;
	float max_length = length.max();
	float resolution = width / max_length;  //每毫米多少个像素


	float *map_img = (float *)malloc(width *height * sizeof(float));
	float *render_img = (float *)malloc(width *height * sizeof(float)); //得到投影渲染后的图像
	mesh_project_xyplane(face_points_idx, all_points, map_img, render_img, resolution, width, height);

	double endTime = clock();
	cout << (endTime - beginTime) << "ms" << endl;

	fstream file;
	file.open("./depth_image.txt", std::ios::out);
	for (int i = 0; i < width *height; i++)
	{
		int y = i / width;
		int x = i % width;
		float z = map_img[i];

		file << x << " " << y << " " << z << std::endl;
		//std::cout << x << " " << y << " " << z << std::endl;
	}
	file.close();

	file.open("./render_image.txt", std::ios::out);
	for (int i = 0; i < width *height; i++)
	{
		int y = i / width;
		int x = i % width;
		float z = render_img[i];

		file << x << " " << y << " " << z << std::endl;
		//std::cout << x << " " << y << " " << z << std::endl;
	}
	file.close();


	free(map_img);
	free(render_img);


	return 0;

}

int mesh_transform_convex_hull(const std::string &file_path, const std::string &save_root)
{

	Triangle_mesh mesh;
	bool flag_io = OpenMesh::IO::read_mesh(mesh, file_path);
	if (false == flag_io)
	{
		std::cout << "load data eror" << std::endl;
		return 0;
	}

	double beginTime = clock();
	int vertic_nums = mesh.n_vertices();
	int face_nums = mesh.n_faces();
	auto points = mesh.points();

	Triangle_mesh::Point min_point(0, 0, 0), max_point(0, 0, 0);
	std::vector<OpenMesh::Vec3f> all_points;
	for (int i = 0; i < vertic_nums; i++)
	{
		min_point = min_point.minimize(points[i]);
		max_point = max_point.maximize(points[i]);
		all_points.push_back(points[i]);
	}

	//得到所有face的点的顶点id号
	std::vector<OpenMesh::Vec3i> face_points_idx;
	get_all_face_points_index(mesh, face_points_idx);

	//减去最小值，使最小点为0 ，然后平移0.5，使最小点为0.5
	transform_points(all_points, min_point);

	Triangle_mesh::Point min_point_ = min_point - min_point;// +Triangle_mesh::Point(0.5, 0.5, 0.5);
	Triangle_mesh::Point max_point_ = max_point - min_point;// +Triangle_mesh::Point(0.5, 0.5, 0.5);


	min_point = Triangle_mesh::Point(1000, 1000, 1000);
	for (int i = 0; i < vertic_nums; i++)
	{
		min_point = min_point.minimize(all_points[i]);
		max_point = max_point.maximize(all_points[i]);
	}


	//平移到中心
	Triangle_mesh::Point trans_point;
	transform_center(all_points, max_point_, min_point_, trans_point);
	trans_point = trans_point - min_point;// +Triangle_mesh::Point(0.5, 0.5, 0.5); //总平移量


	//上面边界平移0.5个距离，那么长度最大值要曾加0.5+0.5 = 1.0个距离，防止边界问题
	//投影到固定尺寸二维平面图像，长宽缩放比一致，这里考虑z轴是垂直空间平面xy轴
	//投影还得考虑中心坐标对齐
	int width = 384;
	int height = 384;
	int depth = 384;
	Triangle_mesh::Point length = max_point_ - min_point_;
	float max_length = length.max();
	float resolution = width / max_length;  //每毫米多少个像素

	int size_ = width * height * depth;
	float *map_img = (float *)malloc(size_ * sizeof(float));
	mesh_voxel_shell_(face_points_idx, all_points, map_img, resolution, width, height, depth);

	double endTime = clock();
	cout << (endTime - beginTime) << "ms" << endl;


	int wh = width * height;
	fstream file;
	file.open("./convex_hull.txt", std::ios::out);
	for (int di = 0; di < size_; di++)
	{
		if (map_img[di] > 0.1)
		{
			int z = di / wh;
			int tmp = di % wh;
			int y = tmp / width;
			int x = tmp % width;

			file << x << " " << y << " " << z << std::endl;
		}

		//std::cout << x << " " << y << " " << z << std::endl;
	}
	file.close();

	free(map_img);


	return 0;

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
		vecFile[i] = "./correct_data.stl";

		mesh_to_image(vecFile[i], save_root);

		mesh_transform_convex_hull(vecFile[i], save_root);

 		std::cout << vecFile[i] << std::endl;
	}

	return 0;
}
