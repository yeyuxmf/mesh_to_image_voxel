#pragma once
#ifndef OPENMESH_WARPPER_H
#define OPENMESH_WARPPER_H
#include"BasicValueTypes.h"
#pragma warning(push)
#pragma warning(disable:4296)
//openmesh
#undef min
#undef max
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>
#pragma warning(pop)

#ifdef USING_X64

typedef double FReal;
namespace OpenMesh {
	struct IntersectTraits
	{
		/// The default coordinate type is OpenMesh::Vec3f.
		typedef Vec3d  Point;

		/// The default normal type is OpenMesh::Vec3f.
		typedef Vec3d  Normal;

		/// The default 1D texture coordinate type is float.
		typedef float  TexCoord1D;
		/// The default 2D texture coordinate type is OpenMesh::Vec2f.
		typedef Vec2f  TexCoord2D;
		/// The default 3D texture coordinate type is OpenMesh::Vec3f.
		typedef Vec3f  TexCoord3D;

		/// The default texture index type
		typedef int TextureIndex;

		/// The default color type is OpenMesh::Vec3uc.
		typedef Vec3uc Color;

#ifndef DOXY_IGNORE_THIS
		VertexTraits{};
		HalfedgeTraits{};
		EdgeTraits{};
		FaceTraits{};
#endif

		VertexAttributes(0);
		HalfedgeAttributes(Attributes::PrevHalfedge);
		EdgeAttributes(0);
		FaceAttributes(0);
	};
}

typedef OpenMesh::TriMesh_ArrayKernelT<OpenMesh::IntersectTraits> Triangle_mesh;
typedef OpenMesh::PolyMesh_ArrayKernelT<OpenMesh::IntersectTraits> Triangle_meshF;
typedef OpenMesh::PolyMesh_ArrayKernelT<OpenMesh::IntersectTraits> Polygon_mesh;
#else
typedef float FReal;
typedef OpenMesh::TriMesh_ArrayKernelT<> Triangle_mesh;
typedef OpenMesh::PolyMesh_ArrayKernelT<> Triangle_meshF;
typedef OpenMesh::PolyMesh_ArrayKernelT<> Polygon_mesh;
#endif
typedef Triangle_mesh::Edge     Edge;
typedef Triangle_mesh::Halfedge Halfedge;
typedef Triangle_mesh::FaceHandle     FaceHandle;
typedef Triangle_mesh::VertexHandle   VertexHandle;
typedef Triangle_mesh::EdgeHandle     EdgeHandle;
typedef Triangle_mesh::HalfedgeHandle HalfedgeHandle;
typedef Triangle_mesh::Point          Point;
typedef Triangle_mesh::Normal         Normal;

#endif  //OPENMESH_WARPPER_H
