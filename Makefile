include makefile.in

INCLUDE = -I$(OPENMESH_INCLUDE_DIR) -Iinclude/ -I$(EIGEN_DIR)
CPPFLAGS = -O3 -fPIC -DEIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS -DEIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET 
LDFLAGS = -O3 -lGL -lGLU
LIB = -lglut -lOpenMeshCored -lOpenMeshToolsd -Wl,-rpath,$(OPENMESH_LIB_DIR)
TARGET = drawMesh
OBJS = objs/main.o objs/curvature.o objs/mesh_features.o objs/image_generation.o objs/decimate.o objs/contour.o objs/highlight.o

default: $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -L$(OPENMESH_LIB_DIR) $(LIB) -o $(TARGET)

objs/main.o: src/main.cpp
	$(CPP) -c $(CPPFLAGS) src/main.cpp -o objs/main.o $(INCLUDE)

objs/curvature.o: src/curvature.cpp
	$(CPP) -c $(CPPFLAGS) src/curvature.cpp -o objs/curvature.o $(INCLUDE)

objs/mesh_features.o: src/mesh_features.cpp
	$(CPP) -c $(CPPFLAGS) src/mesh_features.cpp -o objs/mesh_features.o $(INCLUDE)

objs/image_generation.o: src/image_generation.cpp
	$(CPP) -c $(CPPFLAGS) src/image_generation.cpp -o objs/image_generation.o $(INCLUDE)

objs/decimate.o: src/decimate.cpp
	$(CPP) -c $(CPPFLAGS) src/decimate.cpp -o objs/decimate.o $(INCLUDE)

objs/contour.o: src/contour.cpp
	$(CPP) -c $(CPPFLAGS) src/contour.cpp -o objs/contour.o $(INCLUDE)

objs/highlight.o: src/highlight.cpp
	$(CPP) -c $(CPPFLAGS) src/highlight.cpp -o objs/highlight.o $(INCLUDE)

clean:
	rm $(OBJS) $(TARGET) -f
