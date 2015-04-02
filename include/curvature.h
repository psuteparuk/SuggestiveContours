#ifndef CURVATURE_H
#define CURVATURE_H

#include "mesh_definitions.h"

struct CurvatureInfo {
	OpenMesh::Vec3f directions[2];
	double curvatures[2];
};

void computeCurvature(Mesh &mesh, OpenMesh::VPropHandleT<CurvatureInfo> &curvature);
void computeViewCurvature(Mesh &mesh, OpenMesh::Vec3f camPos, OpenMesh::VPropHandleT<CurvatureInfo> &curvature, OpenMesh::VPropHandleT<double> &viewCurvature, OpenMesh::FPropHandleT<OpenMesh::Vec3f> &viewCurvatureDerivative);

#endif
