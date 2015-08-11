#include "hough_include.h"
#include "hough_detector.h"

int main (int argc, char ** argv)
{
	cimg_usage("Hough circle detector");

	// file input
	const char * file_input = cimg_option("-i", IMAGEPATH "testseq100007.pgm", "Input Image" );
	double sigma = cimg_option("-blur", 1.0, "Standard deviation of gaussian pre-blurring");
	UINT gaussian_window = cimg_option("-window", 5, "Gaussian window size");

	hough_detector detector (file_input, sigma, gaussian_window);
	detector.event_loop();

	return 0;
}
