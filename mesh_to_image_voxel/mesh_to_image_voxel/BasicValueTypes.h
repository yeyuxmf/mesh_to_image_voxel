#pragma once
#ifndef BASIC_VALUE_TYPES_H
#define BASIC_VALUE_TYPES_H




/*! @file BasicValueTypes.h
	@brief 基本的三维数据结构.

	包含基本的数据结构，如三维单精度坐标，双精度坐标，旋转矩阵，包围盒.
*/

// 包含基本的坐标数据结构，如二维双精度PointXY,三维双精度Point
// 包含了stl算法头文件和3*3 4*4旋转矩阵头文件，包围盒头文件，这些文件建议
// 在其他文件中不要单独包含，直接用#include<BasicValueTypes.h>

#pragma warning(push)
//openmesh
#undef min
#undef max
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#define _MATH_DEFINES_DEFINED
#endif
#include<OpenMesh/Core/Geometry/VectorT.hh>
#include<OpenMesh/Core/Mesh/Traits.hh>
#pragma warning(pop)
#include <cstdint>
#include <math.h>
#include <algorithm>

#ifdef USING_X64
typedef double FReal;
#else
typedef float FReal;
#endif

typedef OpenMesh::Vec2d         PointXYd;
typedef OpenMesh::Vec3d         Pointd;
typedef OpenMesh::Vec3d         Normald;

using OpenMesh::VectorT;
#ifndef M_PI_2
#define M_PI_2	1.57079632679489661923
#endif // !


#endif // !BASIC_VALUE_TYPES_H
