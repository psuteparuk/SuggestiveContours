#ifndef HIGHLIGHT_H
#define HIGHLIGHT_H

#include "mesh_definitions.h"

struct HighlightInfo {
  OpenMesh::Vec3f endpoints[2];
  OpenMesh::HalfedgeHandle endedges[2];
  bool hasHighlight;
};

void resetHighlightInfo(Mesh &mesh, OpenMesh::FaceHandle fh, OpenMesh::FPropHandleT<HighlightInfo> &highlight);
void resetAllHighlightInfo(Mesh &mesh, OpenMesh::FPropHandleT<HighlightInfo> &highlight);

#endif
