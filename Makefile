LIBRARIES = -lpthread -lm -lX11 -lstdc++ -fopenmp 
INCLUDES = -I./cimg
CXX_FLAGS = -Wall -W -std=c++11 -Dcimg_use_openmp 

SOURCE = hough_view.cpp hough_detector.cpp main.cpp
OBJECTS = $(SOURCE:.cpp=.o)
EXEC = hough

GCC = gcc

$(EXEC): $(OBJECTS)
	$(GCC) $(CXX_FLAGS) $(INCLUDES) -o $(EXEC) $(OBJECTS) $(LIBRARIES)

.cpp.o:
	$(GCC) $(CXX_FLAGS) $(INCLUDES) -c $< -o $@ $(LIBRARIES)

clean:
	rm -rf $(OBJECTS) $(EXEC)

depend: $(SOURCES)
	makedepend $(INCLUDES) $^
