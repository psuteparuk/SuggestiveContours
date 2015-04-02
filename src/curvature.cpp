#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>
#include "curvature.h"
using namespace OpenMesh;
using namespace Eigen;
using namespace std;

void computeCurvature(Mesh &mesh, OpenMesh::VPropHandleT<CurvatureInfo> &curvature) {

	for (Mesh::VertexIter it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
		// WRITE CODE HERE TO COMPUTE THE CURVATURE AT THE CURRENT VERTEX ----------------------------------------------
                Vec3f pi = mesh.point(it.handle());

		// Compute Normal
                Vec3f normal = mesh.normal(it.handle());
		Vector3d N(normal[0],normal[1],normal[2]);

		// Compute area around v_i
		double total_area_vi = 0.0;
		for (Mesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(it); voh_it; ++voh_it) {
		  total_area_vi += mesh.calc_sector_area(voh_it.handle());
		}

		// Compute Matrix M_i
		Matrix3d M = Matrix3d::Zero();
		for (Mesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(it); voh_it; ++voh_it) {
		  // Compute T_ij
		  Matrix3d IN = Matrix3d::Identity(3,3) - N*N.transpose();
		  Vec3f pj = mesh.point(mesh.to_vertex_handle(voh_it.handle()));
		  Vec3f diffv = pj - pi;
		  Vector3d DiffV(diffv[0],diffv[1],diffv[2]);
		  Vector3d T = IN * DiffV;
		  T.normalize();

		  // Compute k_ij (curvature)
		  double normDiffV = DiffV.norm();
		  double k = 2.0 * N.transpose() * DiffV;
		  k /= (normDiffV * normDiffV);
		  
		  // Compute w_ij
		  double area = 0.0;
		  area += mesh.calc_sector_area(voh_it.handle());
		  area += mesh.calc_sector_area(mesh.opposite_halfedge_handle(voh_it.handle()));
		  double w = area / (2.0*total_area_vi);
		  
		  M += w * k * T * T.transpose();
		}
		
		// Solve eigenvalues
		EigenSolver<Matrix3d> solver(M);
		Matrix3d eigenvectors = solver.pseudoEigenvectors();
		Vector3cd eigenvalues = solver.eigenvalues();
		
		Vector3d v1 = eigenvectors.block(0,0,3,1);
		Vector3d v2 = eigenvectors.block(0,1,3,1);
		Vector3d v3 = eigenvectors.block(0,2,3,1);
		double eig1 = real(eigenvalues(0));
		double eig2 = real(eigenvalues(1));
		double eig3 = real(eigenvalues(2));

		// Find max, min, and 0 eigenvalues
		double m11,m22;
		Vector3d maxv,minv;
		double thresh = 0.01;
		if (v1.cross(N).norm() < thresh) {
		  if (eig2 < eig3) { m11 = eig2; m22 = eig3; maxv = v2; minv = v3; }
		  else { m11 = eig3; m22 = eig2; maxv = v3; minv = v2; }
		}
		else if (v2.cross(N).norm() < thresh) {
		  if (eig3 < eig1) { m11 = eig3; m22 = eig1; maxv = v3; minv = v1; }
		  else { m11 = eig1; m22 = eig3; maxv = v1; minv = v3; }
		}
		else if (v3.cross(N).norm() < thresh) {
		  if (eig1 < eig2) { m11 = eig1; m22 = eig2; maxv = v1; minv = v2; }
		  else { m11 = eig2; m22 = eig1; maxv = v2; minv = v1; }
		}
		maxv.normalize();
		minv.normalize();

		CurvatureInfo info;
		info.curvatures[0] = m22-3*m11;
		info.curvatures[1] = m11-3*m22;
		info.directions[0] = Vec3f(maxv(0),maxv(1),maxv(2));
		info.directions[1] = Vec3f(minv(0),minv(1),minv(2));
		
		mesh.property(curvature,it) = info;
		// -------------------------------------------------------------------------------------------------------------
	}
}

void computeViewCurvature(Mesh &mesh, OpenMesh::Vec3f camPos, OpenMesh::VPropHandleT<CurvatureInfo> &curvature, OpenMesh::VPropHandleT<double> &viewCurvature, OpenMesh::FPropHandleT<OpenMesh::Vec3f> &viewCurvatureDerivative) {
	// WRITE CODE HERE TO COMPUTE CURVATURE IN THE VIEW PROJECTION PROJECTED ON THE TANGENT PLANE ------------------
	// Compute vector to viewer and project onto tangent plane, then // use components in principal directions to find curvature
        for (Mesh::VertexIter it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
	  Vec3f p = mesh.point(it.handle());
	  Vec3f n = mesh.normal(it.handle());
	  Vec3f v = camPos - p;
	  double dist = n | v;
	  Vec3f ortho = n * dist;
	  Vec3f w = v - ortho;
	  w.normalize();

	  CurvatureInfo info = mesh.property(curvature, it);
	  Vec3f T1 = info.directions[0];
	  Vec3f T2 = info.directions[1];
	  double k1 = info.curvatures[0];
	  double k2 = info.curvatures[1];

	  double kw = k1 * (w | T1) * (w | T1) + k2 * (w | T2) * (w | T2);

	  mesh.property(viewCurvature,it) = kw;
	}
	// -------------------------------------------------------------------------------------------------------------

	// We'll use the finite elements piecewise hat method to find per-face gradients of the view curvature
	// CS 348a doesn't cover how to differentiate functions on a mesh (Take CS 468! Spring 2013!) so we provide code here
	
	for (Mesh::FaceIter it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {
		double c[3];
		Vec3f p[3];
		
		Mesh::ConstFaceVertexIter fvIt = mesh.cfv_iter(it);
		for (int i = 0; i < 3; i++) {
			p[i] = mesh.point(fvIt.handle());
			c[i] = mesh.property(viewCurvature,fvIt.handle());
			++fvIt;
		}
		
		Vec3f N = mesh.normal(it.handle());
		double area = mesh.calc_sector_area(mesh.halfedge_handle(it.handle()));

		mesh.property(viewCurvatureDerivative,it) = (N%(p[0]-p[2]))*(c[1]-c[0])/(2*area) + (N%(p[1]-p[0]))*(c[2]-c[0])/(2*area);
	}
}
