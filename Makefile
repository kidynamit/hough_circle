LIBRARIES = -lpthread -lm -lX11 -lstdc++ -fopenmp 
INCLUDES = -I./cimg
CXX_FLAGS = -Wall -W -std=c++11 -Dcimg_use_openmp 

SOURCE = utils.cpp hough_view.cpp kdtree.cpp hough_circles.cpp hough_detector.cpp main.cpp
OBJECTS = $(SOURCE:.cpp=.o)
EXEC = hough

GCC = gcc

$(EXEC): $(OBJECTS)
	$(GCC) $(CXX_FLAGS) $(INCLUDES) -o $(EXEC) $(OBJECTS) $(LIBRARIES)

.cpp.o: softclean
	$(GCC) $(CXX_FLAGS) $(INCLUDES) -c $< -o $@ $(LIBRARIES)

softclean:
	rm -rf $(OBJECTS)

clean: softclean
	rm -rf $(EXEC)

depend: $(SOURCES)
	makedepend $(INCLUDES) $^
