#ifndef  __HOUGH_DETECTOR_H__
#define  __HOUGH_DETECTOR_H__

#include "hough_include.h"

#define DISPLAY_TITLE "Hough Circle Detector"

class hough_detector
{    
private:
	double _sigma;
	int _gaussian_window;
	CImgDisplay _display;
    IMG_LIST_TYPE _images;
	char * _filename;

    UINT _min_radius, _max_radius;
    double _max_threshold, _min_threshold;
    PIXEL_TYPE _hcd_threshold, _hld_threshold;
    void init ();
    void canny_edge_detector ( );
    void hough_circle_detector ( );

    void hough_line_detector (IMG_TYPE & in_image);
	void unmark_line( IMG_TYPE & hough, const int m, const int c);
    void accumulate_line( IMG_TYPE & hough, const int m, const int c);

    int accumulate_circle(IMG_TYPE & hough, const int x, const int y, const int radius );
    void draw_circle(IMG_TYPE & hough, const int x, const int y, const int radius );
	void gaussian_filter ( double **& kernel , double & kernel_sum, const int idx_in = 0, const int idx_out=1);
public:
	hough_detector(const char * file, double sigma=1.0, int gaussian_window=5, double max_thresh=290.0, double min_threshold=200.0);
	~hough_detector();
	void event_loop ();

	static void print_usage ();
};

#endif
