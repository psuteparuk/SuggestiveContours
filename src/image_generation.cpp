#include "image_generation.h"
#include "mesh_features.h"
#include "contour.h"
#include <GL/glut.h>
#include <fstream>
#include <cmath>
#include <set>
#include <map>
#include <list>
using namespace OpenMesh;
using namespace std;

extern VPropHandleT<double> viewCurvature;
extern FPropHandleT<Vec3f> viewCurvatureDerivative;
extern FPropHandleT<ContourInfo> contour;
extern FPropHandleT<bool> chainFlag;

Vec3f toImagePlane(Vec3f point) {
	GLdouble point3DX = point[0], point3DY = point[1], point3DZ = point[2];

	GLdouble modelMatrix[16], projMatrix[16];
	GLint viewport[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
	glGetIntegerv(GL_VIEWPORT, viewport);
	
	GLdouble point2DX, point2DY, point2DZ;
	gluProject(point3DX, point3DY, point3DZ, modelMatrix, projMatrix, viewport, &point2DX, &point2DY, &point2DZ);
	
	return Vec3f(point2DX,point2DY,point2DZ);
}

// Adapted from
// http://stackoverflow.com/questions/1311869/opengl-how-to-determine-if-a-3d-rendered-point-is-occluded-by-other-3d-rende
bool isVisible(Vec3f point) {
	Vec3f projected = toImagePlane(point);

	GLfloat bufDepth = 0.0;
	glReadPixels(static_cast<GLint>( projected[0] ), static_cast<GLint>( projected[1] ),
		     1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &bufDepth);
	std::cout << bufDepth << std::endl;
	GLdouble EPSILON = 0.01;
	return (bufDepth - projected[2]) > -EPSILON; // check sign!
}

void writeImage(Mesh &mesh, int width, int height, string filename, Vec3f camPos, double suggestive_diff_thresh, double suggestive_angle_thresh) {
	ofstream outfile(filename.c_str());
	std::cout << "Start writing to file ..." << std::endl;
	outfile << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
	outfile << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
	outfile << "<svg width=\"5in\" height=\"5in\" viewBox=\"0 0 " << width << ' ' << height << "\" version=\"1.1\" xmlns=\"http:www.w3.org/2000/svg\">\n";
	outfile << "<g stroke=\"black\" fill=\"black\">\n";
	
	// WRITE CODE HERE TO GENERATE A .SVG OF THE MESH --------------------------------------------------------------

	/* Check silhoutte */
	for (Mesh::ConstEdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it)
		if (isFeatureEdge(mesh,*it,camPos)) {
			Mesh::HalfedgeHandle h0 = mesh.halfedge_handle(it,0);
			Mesh::HalfedgeHandle h1 = mesh.halfedge_handle(it,1);
			Vec3f source(mesh.point(mesh.from_vertex_handle(h0)));
			Vec3f target(mesh.point(mesh.from_vertex_handle(h1)));

			//if (!isVisible(source) || !isVisible(target)) continue;
			
			Vec3f p1 = toImagePlane(source);
			Vec3f p2 = toImagePlane(target);
			outfile << "<line ";
			outfile << "x1=\"" << p1[0] << "\" ";
			outfile << "y1=\"" << height-p1[1] << "\" ";
			outfile << "x2=\"" << p2[0] << "\" ";
			outfile << "y2=\"" << height-p2[1] << "\" stroke-width=\"1\" />\n";
		}
	
	/* Render and smooth out suggestive contours */
	ContourInfo info, info1, info2;
	Mesh::HalfedgeHandle heh_0, heh_1;
	Mesh::FaceHandle fh_0, fh_1;
	Mesh::FaceHandle curr, next;
	std::vector< Mesh::FaceHandle > contourLines;
	std::vector<Vec3f> points;

	for (Mesh::ConstFaceIter it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {
	  mesh.property(chainFlag, it) = false;
	}

	for (Mesh::ConstFaceIter it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {
	  info = mesh.property(contour, it);
	  if (!info.hasContour) continue;
	  if (mesh.property(chainFlag, it)) continue; // Checked already

	  heh_0 = info.endedges[0];
	  heh_1 = info.endedges[1];

	  fh_0 = mesh.face_handle(mesh.opposite_halfedge_handle(heh_0));
	  fh_1 = mesh.face_handle(mesh.opposite_halfedge_handle(heh_1));
	  if (mesh.property(contour, fh_0).hasContour && mesh.property(contour, fh_1).hasContour) continue; // Middle of contour lines

	  Vec3f source, target;
	  Vec3f src, mid, end;
	  contourLines.clear();
	  points.clear();
	  curr = it.handle();
	  contourLines.push_back(curr);
	  mesh.property(chainFlag, curr) = true;
	  if (mesh.property(contour, fh_0).hasContour) {
	    next = fh_0;
	    source = info.endpoints[1];
	    target = info.endpoints[0];
	  } else if (mesh.property(contour, fh_1).hasContour) {
	    next = fh_1;
	    source = info.endpoints[0];
	    target = info.endpoints[1];
	  } else {
	    next = fh_0;
	    source = info.endpoints[1];
	    target = info.endpoints[0];
	  }
	  src = toImagePlane(source);
	  end = toImagePlane(target);
	  points.push_back(src);
	  points.push_back(end);

	  while (mesh.property(contour, next).hasContour) {
	    contourLines.push_back(next);
	    info = mesh.property(contour, next);
	    heh_0 = info.endedges[0];
	    heh_1 = info.endedges[1];
	    fh_0 = mesh.face_handle(mesh.opposite_halfedge_handle(heh_0));
	    fh_1 = mesh.face_handle(mesh.opposite_halfedge_handle(heh_1));
	    if (fh_0 == curr) {
	      curr = next;
	      next = fh_1;
	      src = toImagePlane(info.endpoints[1]);
	      points.push_back(src);
	    } else {
	      curr = next;
	      next = fh_0;
	      src = toImagePlane(info.endpoints[0]);
	      points.push_back(src);
	    }

	    mesh.property(chainFlag, curr) = true;
	  }

	  if (points.size() == 2) {
	    src = points[0]; end = points[1];
	    if (abs(src[0]-end[0]) >= 50 || abs(src[1]-end[1]) >= 50) continue;
	    outfile << "<line ";
	    outfile << "x1=\"" << src[0] << "\" ";
	    outfile << "y1=\"" << height-src[1] << "\" ";
	    outfile << "x2=\"" << end[0] << "\" ";
	    outfile << "y2=\"" << height-end[1] << "\" stroke-width=\"1\" />\n";
	    continue;
	  }

	  for (int i = 0; i < points.size()-1; i = i + 2) {
	    src = points[i]; mid = points[i+1]; end = points[i+2];
	    if (abs(src[0]-end[0]) >= 50 || abs(src[1]-end[1]) >= 50) continue;
	    outfile << "<path stroke-width=\"1\" fill=\"none\" d=\"";
	    outfile << "M " << src[0] << "," << height-src[1] << " ";
	    outfile << "Q " << mid[0] << "," << height-mid[1] << " " << end[0] << "," << height-end[1] << "\" />\n";
	  }
	}

	// -------------------------------------------------------------------------------------------------------------
	
	outfile << "</g>\n";
	outfile << "</svg>\n";
	std::cout << "Finished writing to file." << std::endl;

}
