#ifndef MESH_FEATURES_H
#define MESH_FEATURES_H

#include "mesh_definitions.h"

bool isSilhouette(Mesh &mesh, const Mesh::EdgeHandle &e, OpenMesh::Vec3f cameraPos);
bool isSharpEdge(Mesh &mesh, const Mesh::EdgeHandle &e);
bool isFeatureEdge(Mesh &mesh, const Mesh::EdgeHandle &e, OpenMesh::Vec3f cameraPos);

#endif
