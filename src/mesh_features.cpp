#include "mesh_features.h"
using namespace OpenMesh;

bool isSilhouette(Mesh &mesh, const Mesh::EdgeHandle &e, Vec3f cameraPos)  {
    Mesh::HalfedgeHandle heh_0 = mesh.halfedge_handle(e, 0);
    Mesh::HalfedgeHandle heh_1 = mesh.halfedge_handle(e, 1);
    Mesh::FaceHandle fh_0 = mesh.face_handle(heh_0);
    Mesh::FaceHandle fh_1 = mesh.face_handle(heh_1);
    Vec3f n_0 = mesh.normal(fh_0);
    Vec3f n_1 = mesh.normal(fh_1);
    Vec3f pA = mesh.point(mesh.from_vertex_handle(heh_0));
    Vec3f pB = mesh.point(mesh.to_vertex_handle(heh_0));
    Vec3f pC_0, pC_1;
    Mesh::FaceVertexIter fv_it_0 = mesh.fv_iter(fh_0);
    Mesh::FaceVertexIter fv_it_1 = mesh.fv_iter(fh_1);
    Vec3f vt_0, vt_1;
    for (int i = 0; i < 3; ++i) {
      vt_0 = mesh.point(fv_it_0.handle());
      vt_1 = mesh.point(fv_it_1.handle());
      if (vt_0 != pA && vt_0 != pB) pC_0 = vt_0;
      if (vt_1 != pA && vt_1 != pB) pC_1 = vt_1;
      ++fv_it_0;
      ++fv_it_1;
    }
    Vec3f bary_0 = (pA + pB + pC_0) / 3;
    Vec3f bary_1 = (pA + pB + pC_1) / 3;
    Vec3f v_0 = cameraPos - bary_0;
    Vec3f v_1 = cameraPos - bary_1;
    double dot_0 = (n_0 | v_0);
    double dot_1 = (n_1 | v_1);
    return (dot_0 * dot_1 < 0);
}

bool isSharpEdge(Mesh &mesh, const Mesh::EdgeHandle &e) {
    Mesh::HalfedgeHandle heh_0 = mesh.halfedge_handle(e, 0);
    Mesh::HalfedgeHandle heh_1 = mesh.halfedge_handle(e, 1);
    Mesh::FaceHandle fh_0 = mesh.face_handle(heh_0);
    Mesh::FaceHandle fh_1 = mesh.face_handle(heh_1);
    Vec3f n_0 = mesh.normal(fh_0);
    Vec3f n_1 = mesh.normal(fh_1);
    return ((n_0 | n_1) < 0.5);
}

bool isFeatureEdge(Mesh &mesh, const Mesh::EdgeHandle &e, Vec3f cameraPos) {
	return mesh.is_boundary(e) || isSilhouette(mesh,e, cameraPos) || isSharpEdge(mesh,e);
}

