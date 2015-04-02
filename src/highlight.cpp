#include "highlight.h"
using namespace OpenMesh;

void resetHighlightInfo(Mesh &mesh, OpenMesh::FaceHandle fh, OpenMesh::FPropHandleT<HighlightInfo> &highlight) {
  HighlightInfo info;
  info.endpoints[0] = Vec3f();
  info.endpoints[1] = Vec3f();
  info.endedges[0] = HalfedgeHandle();
  info.endedges[1] = HalfedgeHandle();
  info.hasHighlight = false;

  mesh.property(highlight, fh) = info;
}

void resetAllHighlightInfo(Mesh &mesh, OpenMesh::FPropHandleT<HighlightInfo> &highlight) {
  for (Mesh::FaceIter it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {
    resetHighlightInfo(mesh, it.handle(), highlight);
  }
}
