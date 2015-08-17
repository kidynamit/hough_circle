#ifndef  __HOUGH_DETECTOR_H__
#define  __HOUGH_DETECTOR_H__

#include "hough_include.h"
#include "hough_morphology.h"
#include "kdtree.h"

#define DISPLAY_TITLE "Hough Circle Detector"

#define RERENDER { _display.resize ( _display, false).display(_images); }//_display.show ();}

class hough_detector
{    
private:
	double              _sigma;
	int                 _gaussian_window;
	CImgDisplay         _display;
    IMG_LIST_TYPE       _images;
	char *              _filename;
    UINT                _min_radius, _max_radius;
    double              _max_threshold, _min_threshold;
    double              _hcd_threshold, _hcd_low_threshold;
    std::vector <struct kd_node_t> _edge_pixels;
    
	void gaussian_filter ( );
    void init ();
    void canny_edge_detector ( );
    void hough_circle_detector ( );


    double check_circle (const IMG_TYPE & hough, const int x, const int y , const int radius, kdtree & edgetree );
    int accumulate_circle(IMG_TYPE & hough, const int x, const int y, const int radius );
    void draw_circle(IMG_TYPE & hough, const int x, const int y, const int radius);
public:
	hough_detector(const char * file, double sigma=1.0, int gaussian_window=5, double max_thresh=290.0, double min_threshold=200.0);
	~hough_detector();
	
    void event_loop ();
    static void print_usage ();
};

#endif
