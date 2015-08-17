#include "hough_include.h"
#include "hough_detector.h"

int main (int argc, char ** argv)
{
	cimg_usage("Hough circle detector");

	// file input
	const char * file_input = cimg_option("-i", IMAGEPATH "testseq100000.pgm", "Input Image" );
	double sigma = cimg_option("-blur", 2.0, "Standard deviation of gaussian pre-blurring");
	int gaussian_window = cimg_option("-window", 5, "Gaussian window size");
    double max_threshold = cimg_option("-max-intensity", 200.0, "Maximum intensity value. Used for hysteresis");
    double min_threshold = cimg_option("-min-intensity", 100.0, "Minimum intensity value. Used for hysteresis");

	hough_detector detector (file_input, sigma, gaussian_window, max_threshold, min_threshold);
	detector.event_loop();

	return 0;
}
