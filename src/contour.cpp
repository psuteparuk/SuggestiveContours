#include "contour.h"
using namespace OpenMesh;

void resetContourInfo(Mesh &mesh, OpenMesh::FaceHandle fh, OpenMesh::FPropHandleT<ContourInfo> &contour) {
  ContourInfo info;
  info.endpoints[0] = Vec3f();
  info.endpoints[1] = Vec3f();
  info.endedges[0] = HalfedgeHandle();
  info.endedges[1] = HalfedgeHandle();
  info.hasContour = false;

  mesh.property(contour, fh) = info;
}

void resetAllContourInfo(Mesh &mesh, OpenMesh::FPropHandleT<ContourInfo> &contour) {
  for (Mesh::FaceIter it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {
    resetContourInfo(mesh, it.handle(), contour);
  }
}
