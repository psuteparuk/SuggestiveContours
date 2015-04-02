#include <OpenMesh/Core/IO/Options.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <GL/glut.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "curvature.h"
#include "mesh_features.h"
#include "image_generation.h"
#include "decimate.h"
#include "contour.h"
#include "highlight.h"
using namespace std;
using namespace OpenMesh;
using namespace Eigen;

VPropHandleT<double> viewCurvature;
FPropHandleT<Vec3f> viewCurvatureDerivative;
VPropHandleT<CurvatureInfo> curvature;
FPropHandleT<ContourInfo> contour;
FPropHandleT<HighlightInfo> highlight;
FPropHandleT<bool> chainFlag;
Mesh mesh;

bool leftDown = false, rightDown = false, middleDown = false;
int lastPos[2];
float cameraPos[4] = {0,0,4,1};
Vec3f up, pan;
int windowWidth = 640, windowHeight = 480;
float axisAngle = 0.0, animateAngle = 0.0;
bool showSurface = true, showAxes = true, showCurvature = false;
bool showNormals = false, showMesh = false, showLighting = true;
bool showSuggestiveContours = true, showSuggestiveHighlights = false;
bool filterLines = false;
bool showRedGreen = false;

float specular[] = { 1.0, 1.0, 1.0, 1.0 };
float shininess[] = { 50.0 };

#define PI 3.14159265

double mark_length = 0.02;
double suggestive_diff_thresh_min = 20.0;
double suggestive_diff_thresh_max = 600.0;
double suggestive_diff_thresh_incr = 20.0;
double suggestive_angle_thresh_min = 6.0; // degree
double suggestive_angle_thresh_max = 20.0;
double suggestive_angle_thresh_incr = 2.0;
unsigned int contour_length_thresh_min = 1;
unsigned int contour_length_thresh_max = 5;
unsigned int contour_length_thresh_incr = 1;

double suggestive_diff_thresh = suggestive_diff_thresh_min;
double suggestive_angle_thresh = suggestive_angle_thresh_min;
unsigned int contour_length_thresh = contour_length_thresh_min;

void clearChainFlag() {
        for (Mesh::ConstFaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
	  mesh.property(chainFlag, f_it) = false;
	}
}

void filterContourLines() {
	ContourInfo info;
	Mesh::HalfedgeHandle heh_0, heh_1;
	Mesh::FaceHandle fh_0, fh_1;
	Mesh::FaceHandle curr, next;
	std::vector< Mesh::FaceHandle > contourLines;
	
	clearChainFlag();

	for (Mesh::ConstFaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
	  info = mesh.property(contour, f_it);
	  if (!info.hasContour) continue;
	  if (mesh.property(chainFlag, f_it)) continue; // Checked already

	  heh_0 = info.endedges[0];
	  heh_1 = info.endedges[1];
	  
	  fh_0 = mesh.face_handle(mesh.opposite_halfedge_handle(heh_0));
	  fh_1 = mesh.face_handle(mesh.opposite_halfedge_handle(heh_1));
	  if (mesh.property(contour, fh_0).hasContour && mesh.property(contour, fh_1).hasContour) continue; // Middle of contour lines
	  
	  contourLines.clear();
	  curr = f_it.handle();
	  contourLines.push_back(curr);
	  mesh.property(chainFlag, curr) = true;
	  if (mesh.property(contour, fh_0).hasContour) next = fh_0;
	  else if (mesh.property(contour, fh_1).hasContour) next = fh_1;
	  else next = fh_0;

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
	    } else {
	      curr = next;
	      next = fh_0;
	    }
	    mesh.property(chainFlag, curr) = true;
	  }

	  if (contourLines.size() <= contour_length_thresh) {
	    for (int i = 0; i < contourLines.size(); ++i) {
	      resetContourInfo(mesh, contourLines[i], contour);
	    }
	  }
	}
}

void filterHighlightLines() {
	HighlightInfo info;
	Mesh::HalfedgeHandle heh_0, heh_1;
	Mesh::FaceHandle fh_0, fh_1;
	Mesh::FaceHandle curr, next;
	std::vector< Mesh::FaceHandle > highlightLines;
	
	clearChainFlag();

	for (Mesh::ConstFaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
	  info = mesh.property(highlight, f_it);
	  if (!info.hasHighlight) continue;
	  if (mesh.property(chainFlag, f_it)) continue; // Checked already

	  heh_0 = info.endedges[0];
	  heh_1 = info.endedges[1];
	  
	  fh_0 = mesh.face_handle(mesh.opposite_halfedge_handle(heh_0));
	  fh_1 = mesh.face_handle(mesh.opposite_halfedge_handle(heh_1));
	  if (mesh.property(highlight, fh_0).hasHighlight && mesh.property(highlight, fh_1).hasHighlight) continue; // Middle of highlight lines
	  
	  highlightLines.clear();
	  curr = f_it.handle();
	  highlightLines.push_back(curr);
	  mesh.property(chainFlag, curr) = true;
	  if (mesh.property(highlight, fh_0).hasHighlight) next = fh_0;
	  else if (mesh.property(highlight, fh_1).hasHighlight) next = fh_1;
	  else next = fh_0;

	  while (mesh.property(highlight, next).hasHighlight) {
	    highlightLines.push_back(next);
	    info = mesh.property(highlight, next);
	    heh_0 = info.endedges[0];
	    heh_1 = info.endedges[1];
	    fh_0 = mesh.face_handle(mesh.opposite_halfedge_handle(heh_0));
	    fh_1 = mesh.face_handle(mesh.opposite_halfedge_handle(heh_1));
	    if (fh_0 == curr) {
	      curr = next;
	      next = fh_1;
	    } else {
	      curr = next;
	      next = fh_0;
	    }
	    mesh.property(chainFlag, curr) = true;
	  }

	  if (highlightLines.size() <= contour_length_thresh) {
	    for (int i = 0; i < highlightLines.size(); ++i) {
	      resetHighlightInfo(mesh, highlightLines[i], highlight);
	    }
	  }
	}
}

void renderSuggestiveContours(Vec3f actualCamPos) { // use this camera position to account for panning etc.
	glColor3f(0,0,0);
	if (showRedGreen) glColor3f(1,0,0);
	
	resetAllContourInfo(mesh, contour);

	// RENDER SUGGESTIVE CONTOURS HERE -----------------------------------------------------------------------------
	for (Mesh::ConstFaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
	  Vec3f Dw = mesh.property(viewCurvatureDerivative, f_it);

	  Vec3f n = mesh.normal(f_it.handle());
	  Vec3f pA, pB, pC;
	  Mesh::VertexHandle pAh, pBh, pCh;
	  Mesh::HalfedgeHandle heh_1, heh_2, heh_3;
	  Mesh::ConstFaceHalfedgeIter cfh_it = mesh.cfh_iter(f_it.handle());
	  heh_1 = cfh_it.handle();
	  heh_2 = (++cfh_it).handle();
	  heh_3 = (++cfh_it).handle();
	  pAh = mesh.from_vertex_handle(heh_1);
	  pBh = mesh.from_vertex_handle(heh_2);
	  pCh = mesh.from_vertex_handle(heh_3);
	  pA = mesh.point(pAh);
	  pB = mesh.point(pBh);
	  pC = mesh.point(pCh);

	  double kwA = mesh.property(viewCurvature, pAh);
	  double kwB = mesh.property(viewCurvature, pBh);
	  double kwC = mesh.property(viewCurvature, pCh);
	  
	  Vec3f p1, p2;
	  Mesh::HalfedgeHandle h1, h2;
	  double t;
	  bool hasone = false, hastwo = false;
	  if (kwA * kwB < 0) {
	    t = kwA / (kwA - kwB);
	    p1 = pA + (pB - pA) * t;
	    h1 = heh_1;
	    hasone = true;
	  }
	  if (kwB * kwC < 0) {
	    t = kwB / (kwB - kwC);
	    if (!hasone) {
	      p1 = pB + (pC - pB) * t;
	      h1 = heh_2;
	    } else {
	      p2 = pB + (pC - pB) * t;
	      h2 = heh_2;
	      hastwo = true;
	    }
	    hasone = true;
	  }
	  if (kwC * kwA < 0) {
	    t = kwC / (kwC - kwA);
	    p2 = pC + (pA - pC) * t;
	    h2 = heh_3;
	    hastwo = true;
	  }
	  
	  if (!hastwo) continue;

	  Vec3f bary = (pA + pB + pC) / 3;
	  Vec3f v = actualCamPos - bary;
	  double cos_theta = (n | v) / v.norm();
	  if (cos_theta > cos(suggestive_angle_thresh*PI/180)) continue;

	  Vec3f w1 = actualCamPos - p1; w1.normalize();
	  Vec3f w2 = actualCamPos - p2; w2.normalize();
	  double d1 = (Dw | w1);
	  double d2 = (Dw | w2);
	  if (d1 <= suggestive_diff_thresh && d2 <= suggestive_diff_thresh) continue;
	  if (d1 > suggestive_diff_thresh && d2 <= suggestive_diff_thresh) {
	    t = (d1 - suggestive_diff_thresh) / (d1 - d2);
	    p2 = p1 + (p2 - p1) * t;
	  } else if (d1 <= suggestive_diff_thresh && d2 > suggestive_diff_thresh) {
	    t = (d2 - suggestive_diff_thresh) / (d2 - d1);
	    p1 = p2 + (p1 - p2) * t;
	  }

	  ContourInfo info;
	  info.endpoints[0] = p1;
	  info.endpoints[1] = p2;
	  info.endedges[0] = h1;
	  info.endedges[1] = h2;
	  info.hasContour = true;
	  
	  mesh.property(contour, f_it) = info;
	}

	if (filterLines) filterContourLines();

	glBegin(GL_LINES);
	for (Mesh::ConstFaceIter it = mesh.faces_begin(); it !=	mesh.faces_end(); ++it) {
	  ContourInfo info = mesh.property(contour, it);
	  if (info.hasContour) {
	    Vec3f p1 = info.endpoints[0];
	    Vec3f p2 = info.endpoints[1];
	    glVertex3f(p1[0],p1[1],p1[2]);
	    glVertex3f(p2[0],p2[1],p2[2]);
	  }
	}
	glEnd();
	// -------------------------------------------------------------------------------------------------------------
}

void renderSuggestiveHighlights(Vec3f actualCamPos) { // use this camera position to account for panning etc.
	if (showSurface) glColor3f(1,1,1);
	else glColor3f(0.5,0.5,0.5);
	if (showRedGreen) glColor3f(0,1,0);
	
	resetAllHighlightInfo(mesh, highlight);

	// RENDER SUGGESTIVE HIGHLIGHTS HERE -----------------------------------------------------------------------------
	for (Mesh::ConstFaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
	  Vec3f Dw = mesh.property(viewCurvatureDerivative, f_it);

	  Vec3f n = mesh.normal(f_it.handle());
	  Vec3f pA, pB, pC;
	  Mesh::VertexHandle pAh, pBh, pCh;
	  Mesh::HalfedgeHandle heh_1, heh_2, heh_3;
	  Mesh::ConstFaceHalfedgeIter cfh_it = mesh.cfh_iter(f_it.handle());
	  heh_1 = cfh_it.handle();
	  heh_2 = (++cfh_it).handle();
	  heh_3 = (++cfh_it).handle();
	  pAh = mesh.from_vertex_handle(heh_1);
	  pBh = mesh.from_vertex_handle(heh_2);
	  pCh = mesh.from_vertex_handle(heh_3);
	  pA = mesh.point(pAh);
	  pB = mesh.point(pBh);
	  pC = mesh.point(pCh);

	  double kwA = mesh.property(viewCurvature, pAh);
	  double kwB = mesh.property(viewCurvature, pBh);
	  double kwC = mesh.property(viewCurvature, pCh);
	  
	  Vec3f p1, p2;
	  Mesh::HalfedgeHandle h1, h2;
	  double t;
	  bool hasone = false, hastwo = false;
	  if (kwA * kwB < 0) {
	    t = kwA / (kwA - kwB);
	    p1 = pA + (pB - pA) * t;
	    h1 = heh_1;
	    hasone = true;
	  }
	  if (kwB * kwC < 0) {
	    t = kwB / (kwB - kwC);
	    if (!hasone) {
	      p1 = pB + (pC - pB) * t;
	      h1 = heh_2;
	    } else {
	      p2 = pB + (pC - pB) * t;
	      h2 = heh_2;
	      hastwo = true;
	    }
	    hasone = true;
	  }
	  if (kwC * kwA < 0) {
	    t = kwC / (kwC - kwA);
	    p2 = pC + (pA - pC) * t;
	    h2 = heh_3;
	    hastwo = true;
	  }
	  
	  if (!hastwo) continue;

	  Vec3f bary = (pA + pB + pC) / 3;
	  Vec3f v = actualCamPos - bary;
	  double cos_theta = (n | v) / v.norm();
	  if (cos_theta > cos(suggestive_angle_thresh*PI/180)) continue;

	  Vec3f w1 = actualCamPos - p1; w1.normalize();
	  Vec3f w2 = actualCamPos - p2; w2.normalize();
	  double d1 = (Dw | w1);
	  double d2 = (Dw | w2);
	  if (d1 >= -suggestive_diff_thresh*2 && d2 >= -suggestive_diff_thresh*2) continue;
	  if (d1 < -suggestive_diff_thresh*2 && d2 >= -suggestive_diff_thresh*2) {
	    t = (d1 + suggestive_diff_thresh*2) / (d1 - d2);
	    p2 = p1 + (p2 - p1) * t;
	  } else if (d1 >= -suggestive_diff_thresh*2 && d2 < -suggestive_diff_thresh*2) {
	    t = (d2 + suggestive_diff_thresh*2) / (d2 - d1);
	    p1 = p2 + (p1 - p2) * t;
	  }

	  HighlightInfo info;
	  info.endpoints[0] = p1;
	  info.endpoints[1] = p2;
	  info.endedges[0] = h1;
	  info.endedges[1] = h2;
	  info.hasHighlight = true;
	  
	  mesh.property(highlight, f_it) = info;
	}

	if (filterLines) filterHighlightLines();

	glBegin(GL_LINES);
	for (Mesh::ConstFaceIter it = mesh.faces_begin(); it !=	mesh.faces_end(); ++it) {
	  HighlightInfo info = mesh.property(highlight, it);
	  if (info.hasHighlight) {
	    Vec3f p1 = info.endpoints[0];
	    Vec3f p2 = info.endpoints[1];
	    glVertex3f(p1[0],p1[1],p1[2]);
	    glVertex3f(p2[0],p2[1],p2[2]);
	  }
	}
	glEnd();
	// -------------------------------------------------------------------------------------------------------------
}

void renderMesh() {
	if (!showSurface) glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE); // render regardless to remove hidden lines
	
	glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0, GL_POSITION, cameraPos);
	if (!showLighting) glDisable(GL_LIGHTING);

	glDepthRange(0.001,1);
	glEnable(GL_NORMALIZE);
 
	// WRITE CODE HERE TO RENDER THE TRIANGLES OF THE MESH ---------------------------------------------------------
	for (Mesh::ConstFaceIter it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {
	  Vec3f nA, nB, nC;
	  Vec3f pA, pB, pC;
	  Mesh::ConstFaceVertexIter cfv_it = mesh.cfv_iter(it.handle());
	  nA = mesh.normal(cfv_it.handle());
	  pA = mesh.point(cfv_it.handle());
	  nB = mesh.normal((++cfv_it).handle());
	  pB = mesh.point(cfv_it.handle());
	  nC = mesh.normal((++cfv_it).handle());
	  pC = mesh.point(cfv_it.handle());
	  glBegin(GL_TRIANGLES);
	  glColor3f(0.7,0.7,0.7);
	  glNormal3f(nA[0],nA[1],nA[2]);
	  glVertex3f(pA[0],pA[1],pA[2]);
	  glNormal3f(nB[0],nB[1],nB[2]);
	  glVertex3f(pB[0],pB[1],pB[2]);
	  glNormal3f(nC[0],nC[1],nC[2]);
	  glVertex3f(pC[0],pC[1],pC[2]);
	  glEnd();
	}
	// -------------------------------------------------------------------------------------------------------------
	
	if (!showSurface) glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
	
	glDisable(GL_LIGHTING);
	glDepthRange(0,0.999);
	
	Vec3f actualCamPos(cameraPos[0]+pan[0],cameraPos[1]+pan[1],cameraPos[2]+pan[2]);
	if (showSuggestiveContours) renderSuggestiveContours(actualCamPos);
	if (showSuggestiveHighlights) renderSuggestiveHighlights(actualCamPos);
	
	// We'll be nice and provide you with code to render feature edges below
	glBegin(GL_LINES);
	glColor3f(0,0,0);
	glLineWidth(2.0f);
	for (Mesh::ConstEdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it)
		if (showMesh || isFeatureEdge(mesh,*it,actualCamPos)) {
			Mesh::HalfedgeHandle h0 = mesh.halfedge_handle(it,0);
			Mesh::HalfedgeHandle h1 = mesh.halfedge_handle(it,1);
			Vec3f source(mesh.point(mesh.from_vertex_handle(h0)));
			Vec3f target(mesh.point(mesh.from_vertex_handle(h1)));
			glVertex3f(source[0],source[1],source[2]);
			glVertex3f(target[0],target[1],target[2]);
		}
	glEnd();
	
	if (showCurvature) {
		// WRITE CODE HERE TO RENDER THE PRINCIPAL DIRECTIONS YOU COMPUTED ---------------------------------------------
	  glBegin(GL_LINES);
	  glLineWidth(1);
	  for (Mesh::ConstVertexIter it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
	    CurvatureInfo info = mesh.property(curvature,it);
	    Vec3f maxdir = info.directions[0];
	    Vec3f mindir = info.directions[1];
	    Vec3f p = mesh.point(it.handle());
	    Vec3f maxd1 = p - maxdir*mark_length;
	    Vec3f maxd2 = p + maxdir*mark_length;
	    Vec3f mind1 = p - mindir*mark_length;
	    Vec3f mind2 = p + mindir*mark_length;
	    glColor3f(0,0,1); glVertex3f(maxd1[0],maxd1[1],maxd1[2]); glVertex3f(maxd2[0], maxd2[1], maxd2[2]);
	    glColor3f(1,0,0); glVertex3f(mind1[0],mind1[1],mind1[2]); glVertex3f(mind2[0], mind2[1], mind2[2]);
	  }
	  glEnd();
		// -------------------------------------------------------------------------------------------------------------
	}
	
	if (showNormals) {
		glBegin(GL_LINES);
		glColor3f(0,1,0);
		for (Mesh::ConstVertexIter it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
			Vec3f n = mesh.normal(it.handle());
			Vec3f p = mesh.point(it.handle());
			Vec3f d = p + n*2*mark_length;
			glVertex3f(p[0],p[1],p[2]);
			glVertex3f(d[0],d[1],d[2]);
		}
		glEnd();
	}
	
	glDepthRange(0,1);
}

void display() {
	glClearColor(1,1,1,1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, shininess);
	glEnable(GL_LIGHT0);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0,0,windowWidth,windowHeight);
	
	float ratio = (float)windowWidth / (float)windowHeight;
	gluPerspective(50, ratio, 1, 1000); // 50 degree vertical viewing angle, zNear = 1, zFar = 1000
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cameraPos[0]+pan[0], cameraPos[1]+pan[1], cameraPos[2]+pan[2], pan[0], pan[1], pan[2], up[0], up[1], up[2]);
	
	// Draw mesh
	renderMesh();

	// Draw axes
	if (showAxes) {
		glDisable(GL_LIGHTING);
		glBegin(GL_LINES);
		glLineWidth(1);
			glColor3f(1,0,0); glVertex3f(0,0,0); glVertex3f(1,0,0); // x axis
			glColor3f(0,1,0); glVertex3f(0,0,0); glVertex3f(0,1,0); // y axis
			glColor3f(0,0,1); glVertex3f(0,0,0); glVertex3f(0,0,1); // z axis
		glEnd(/*GL_LINES*/);
	}

	glutSwapBuffers();
}

void mouse(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON) leftDown = (state == GLUT_DOWN);
	else if (button == GLUT_RIGHT_BUTTON) rightDown = (state == GLUT_DOWN);
	else if (button == GLUT_MIDDLE_BUTTON) middleDown = (state == GLUT_DOWN);
	
	lastPos[0] = x;
	lastPos[1] = y;
}

void mouseMoved(int x, int y) {
	int dx = x - lastPos[0];
	int dy = y - lastPos[1];
	Vec3f curCamera(cameraPos[0],cameraPos[1],cameraPos[2]);
	Vec3f curCameraNormalized = curCamera.normalized();
	Vec3f right = up % curCameraNormalized;

	if (leftDown) {
		// Assume here that up vector is (0,1,0)
		Vec3f newPos = curCamera - 2*(float)((float)dx/(float)windowWidth) * right + 2*(float)((float)dy/(float)windowHeight) * up;
		newPos = newPos.normalized() * curCamera.length();
		
		up = up - (up | newPos) * newPos / newPos.sqrnorm();
		up.normalize();
		
		for (int i = 0; i < 3; i++) cameraPos[i] = newPos[i];
	}
	else if (rightDown) for (int i = 0; i < 3; i++) cameraPos[i] *= pow(1.1,dy*.1);
	else if (middleDown) {
		pan += -2*(float)((float)dx/(float)windowWidth) * right + 2*(float)((float)dy/(float)windowHeight) * up;
	}

	
	lastPos[0] = x;
	lastPos[1] = y;
	
	Vec3f actualCamPos(cameraPos[0]+pan[0],cameraPos[1]+pan[1],cameraPos[2]+pan[2]);
	computeViewCurvature(mesh,actualCamPos,curvature,viewCurvature,viewCurvatureDerivative);
	
	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
	Vec3f actualCamPos(cameraPos[0]+pan[0],cameraPos[1]+pan[1],cameraPos[2]+pan[2]);

	if (showSuggestiveContours || showSuggestiveHighlights) {
	  if (key == '+' && suggestive_angle_thresh < suggestive_angle_thresh_max) suggestive_angle_thresh += suggestive_angle_thresh_incr;
	  else if (key == '+' && suggestive_angle_thresh == suggestive_angle_thresh_max) std::cout << "Max threshold" << std::endl;
	  else if (key == '-' && suggestive_angle_thresh > suggestive_angle_thresh_min) suggestive_angle_thresh -= suggestive_angle_thresh_incr;
	  else if (key == '-' && suggestive_angle_thresh == suggestive_angle_thresh_min) std::cout << "Min threshold" << std::endl;
	}

	if (key == 's' || key == 'S') showSurface = !showSurface;
	else if (key == 'a' || key == 'A') showAxes = !showAxes;
	else if (key == 'c' || key == 'C') showCurvature = !showCurvature;
	else if (key == 'n' || key == 'N') showNormals = !showNormals;
	else if (key == 'm' || key == 'M') showMesh = !showMesh;
	else if (key == 'l' || key == 'L') showLighting = !showLighting;
	else if (key == 'g' || key == 'G') showSuggestiveContours = !showSuggestiveContours;
	else if (key == 'h' || key == 'H') showSuggestiveHighlights = !showSuggestiveHighlights;
	else if (key == 'f' || key == 'F') filterLines = !filterLines;
	else if (key == ' ') showRedGreen = !showRedGreen;
	else if (key == 'w' || key == 'W') writeImage(mesh, windowWidth, windowHeight, "renderedImage.svg", actualCamPos, suggestive_diff_thresh, suggestive_angle_thresh);
	else if (key == 'q' || key == 'Q') exit(0);

	glutPostRedisplay();
}

void keyspecial(int key, int x, int y) {
        if (showSuggestiveContours || showSuggestiveHighlights) {
	  if (key == GLUT_KEY_UP && suggestive_diff_thresh < suggestive_diff_thresh_max) suggestive_diff_thresh += suggestive_diff_thresh_incr;
	  else if (key == GLUT_KEY_UP && suggestive_diff_thresh == suggestive_diff_thresh_max) std::cout << "Max threshold" << std::endl;
	  else if (key == GLUT_KEY_DOWN && suggestive_diff_thresh > suggestive_diff_thresh_min) suggestive_diff_thresh -= suggestive_diff_thresh_incr;
	  else if (key == GLUT_KEY_DOWN && suggestive_diff_thresh == suggestive_diff_thresh_min) std::cout << "Min threshold" << std::endl;
	}
	
	if (filterLines) {
	  if (key == GLUT_KEY_PAGE_UP && contour_length_thresh < contour_length_thresh_max) contour_length_thresh += contour_length_thresh_incr;
	  else if (key == GLUT_KEY_PAGE_UP && contour_length_thresh == contour_length_thresh_max) std::cout << "Max threshold" << std::endl;
	  else if (key == GLUT_KEY_PAGE_DOWN && contour_length_thresh > contour_length_thresh_min) contour_length_thresh -= contour_length_thresh_incr;
	  else if (key == GLUT_KEY_PAGE_DOWN && contour_length_thresh == contour_length_thresh_min) std::cout << "Min threshold" << std::endl;
	}

	glutPostRedisplay();
}

void reshape(int width, int height) {
	windowWidth = width;
	windowHeight = height;
	glutPostRedisplay();
}

int main(int argc, char** argv) {
	if (argc < 2) {
		cout << "Usage: " << argv[0] << " mesh_filename\n";
		exit(0);
	}
	
	IO::Options opt;
	opt += IO::Options::VertexNormal;
	opt += IO::Options::FaceNormal;
	
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	
	cout << "Reading from file " << argv[1] << "...\n";
	if ( !IO::read_mesh(mesh, argv[1], opt )) {
		cout << "Read failed.\n";
		exit(0);
	}

	cout << "Mesh stats:\n";
	cout << '\t' << mesh.n_vertices() << " vertices.\n";
	cout << '\t' << mesh.n_edges() << " edges.\n";
	cout << '\t' << mesh.n_faces() << " faces.\n";
	
	simplify(mesh,.5f);
	
	mesh.update_normals();
	
	mesh.add_property(viewCurvature);
	mesh.add_property(viewCurvatureDerivative);
	mesh.add_property(curvature);
	mesh.add_property(contour);
	mesh.add_property(highlight);
	mesh.add_property(chainFlag);
	
	// Move center of mass to origin
	Vec3f center(0,0,0);
	for (Mesh::ConstVertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) center += mesh.point(vIt);
	center /= mesh.n_vertices();
	for (Mesh::VertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) mesh.point(vIt) -= center;

	// Fit in the unit sphere
	float maxLength = 0;
	for (Mesh::ConstVertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) maxLength = max(maxLength, mesh.point(vIt).length());
	for (Mesh::VertexIter vIt = mesh.vertices_begin(); vIt != mesh.vertices_end(); ++vIt) mesh.point(vIt) /= maxLength;
	
	computeCurvature(mesh,curvature);

	up = Vec3f(0,1,0);
	pan = Vec3f(0,0,0);
	
	Vec3f actualCamPos(cameraPos[0]+pan[0],cameraPos[1]+pan[1],cameraPos[2]+pan[2]);
	computeViewCurvature(mesh,actualCamPos,curvature,viewCurvature,viewCurvatureDerivative);

	glutInit(&argc, argv); 
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); 
	glutInitWindowSize(windowWidth, windowHeight); 
	glutCreateWindow(argv[0]);

	glutDisplayFunc(display);
	glutMotionFunc(mouseMoved);
	glutMouseFunc(mouse);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(keyspecial);

	glutMainLoop();
	
	return 0;
}
