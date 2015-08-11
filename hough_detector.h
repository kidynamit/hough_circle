#ifndef  __HOUGH_DETECTOR_H__
#define  __HOUGH_DETECTOR_H__

#include "hough_include.h"

class hough_detector
{    
private:
	double _sigma;
	UINT _gaussian_window;
	CImgDisplay _display, _gaussian_display;
	IMG_TYPE _image;
	char * _filename;
	void init ();
	void gaussian_filter ();
public:
	hough_detector(const char * file, double sigma=1.0, UINT gaussian_window=5);
	~hough_detector();
	void event_loop ();

	static void print_usage ();
};

#endif 

