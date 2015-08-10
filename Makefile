LIBRARIES = -lpthread -lm -lX11
INCLUDES = -I./cimg

hough:
	g++ -o hough main.cpp $(INCLUDES) $(LIBRARIES)

clean:
	rm -rf main.o hough
