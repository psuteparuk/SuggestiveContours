#ifndef CONTOUR_H
#define CONTOUR_H

#include "mesh_definitions.h"

struct ContourInfo {
  OpenMesh::Vec3f endpoints[2];
  OpenMesh::HalfedgeHandle endedges[2];
  bool hasContour;
};

void resetContourInfo(Mesh &mesh, OpenMesh::FaceHandle fh, OpenMesh::FPropHandleT<ContourInfo> &contour);
void resetAllContourInfo(Mesh &mesh, OpenMesh::FPropHandleT<ContourInfo> &contour);

#endif
