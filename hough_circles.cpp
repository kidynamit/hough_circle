#include "hough_detector.h"

void hough_detector::hough_circle_detector ()
{
	if ( _edge_pixels.empty() )
	{
		LOG(ERROR, "no edge pixels");
		return;
	}

	struct kd_node_t * edge_pixels_arr = new struct kd_node_t [ _edge_pixels.size () ];
	for (UINT i = 0; i < _edge_pixels.size() ; i++ )
	{
		edge_pixels_arr[i] = _edge_pixels[i];
	}

	// construct kdtree
	kdtree edge_tree (edge_pixels_arr, 2);

	std::stringstream ss, debug_info;
	ss << "Applying Hough Circle Detection: \t\t" << _hcd_threshold <<" - "  <<  _min_radius<< " , "<< _max_radius; 
	LOG(INFO, ss.str().c_str() );

	_display.set_title ( (DISPLAY_TITLE " ... Applying Hough Circle Detection") );
	
	// _display.close();
	const int width = _images[0].width();
	const int height = _images[0].height();
	
	// construct list 
	IMG_TYPE accumulator (width, height );
	cimg_forXY(accumulator, x, y) 
    {
        accumulator(x, y) = NULL_PIXEL;
    }

	IMG_TYPE hough (width, height);
	for ( int radius = _min_radius; radius < (int)_max_radius ; radius ++ ) 
	{
		
		cimg_forXY(hough, x, y) hough(x, y) = NULL_PIXEL;
		int cum_circumference = 0;
		int total = 0;
		for ( int x = 0 ;  x < width ; x ++ ) 
		{
			for ( int y = 0; y < height ; y++ )
			{
				if ( _images[3](x, y) == STRONG_PIXEL)
				{            
					cum_circumference += accumulate_circle ( hough, x, y, radius);
					total ++;
				}
			}
		}

        cimg_forXY(accumulator, x, y) if ( hough(x, y) > accumulator (x, y) ) accumulator(x, y) = hough(x, y);

        double circumference = cum_circumference  / (double) total;

		double high_threshold ; 
        double low_threshold ;
		// double threshold = _hcd_threshold ;
		for ( int x = 0;  x < width ; x++ ) 
		{
			for ( int y = 0; y < height ; y++ )
			{
                int bound_x = MIN(x + radius, width - 1) - MAX(x - radius, 0);
                int bound_y = MIN(y + radius, height - 1) - MAX(y - radius, 0);

                double occlusion_factor = ((double) 4.0 * radius * radius) / ((double)( bound_x * bound_y )) ;
                high_threshold = circumference / (_hcd_threshold * occlusion_factor) ;
                low_threshold = (_hcd_low_threshold * circumference) / (_hcd_threshold * occlusion_factor) ;
                if ( ( hough(x, y)) > high_threshold )
				{
                    draw_circle ( accumulator, x, y, radius);
				}
                else if ( hough(x, y) > low_threshold )
                {
                    // TODO check for the area difference of the circle drawn with the edge pixels present
                    if ( check_circle ( accumulator, x, y, radius, edge_tree ) > 0.0003 * low_threshold  )
                    {
                       draw_circle ( accumulator, x, y, radius);
                    }
                }
			}
		}
	}
    
    delete [] edge_pixels_arr;

	_images[4] = accumulator;
	
	_display.set_title ( (DISPLAY_TITLE " ... Ready!") );
	cimg::wait (1000);
	RERENDER;
	_display.set_title ( (DISPLAY_TITLE) );
}

int hough_detector::accumulate_circle (IMG_TYPE & hough, const int x, const int y , const int radius )
{
	const int width = hough.width ();
	const int height = hough.height ();

	int r = radius;
	int c = 0;
	int num_pixels_set = 0;
	int decider_2 = 1 - r ;
	while (r >= c)
	{
		for ( int i = -1; i <= 1 ; i += 2)
		{
			for ( int j = -1; j <= 1; j += 2) 
			{
				int idx_x = x + r * i;
				int idx_y = y + c * j;
				
				if ( idx_x >= 0 && idx_x < width && idx_y >= 0 && idx_y < height ) 
				{	
					hough(idx_x, idx_y) = hough(idx_x, idx_y) + 1;
					num_pixels_set ++;
				}
				idx_x = x + c * i;
				idx_y = y + r * j;
				
				if ( idx_x >= 0 && idx_x < width && idx_y >= 0 && idx_y < height ) 
				{
					hough(idx_x, idx_y) = hough(idx_x, idx_y) + 1;
					num_pixels_set ++;
				}
			}
		}

		c++; 
		if ( decider_2 <= 0 )
		{
			decider_2 += 2 * c + 1;
		}
		else
		{
			r-- ;
			decider_2 += 2 * ( c - r ) + 1;
		}

	}
	return num_pixels_set; 
}

void hough_detector::draw_circle (IMG_TYPE & hough, const int x, const int y , const int radius)
{
	const int width = hough.width ();
	const int height = hough.height ();
    
	int r = radius;
	int c = 0;
	int decider_2 = 1 - r ;
	while (r >= c)
	{
		for ( int i = -1; i <= 1 ; i += 2)
		{
			for ( int j = -1; j <= 1; j += 2) 
			{
				int idx_x = x + r * i;
				int idx_y = y + c * j;
				
				if ( idx_x >= 0 && idx_x < width && idx_y >= 0 && idx_y < height ) 
				{
					hough(idx_x, idx_y) = STRONG_PIXEL;
				}
				idx_x = x + c * i;
				idx_y = y + r * j;
				
				if ( idx_x >= 0 && idx_x < width && idx_y >= 0 && idx_y < height ) 
				{
					hough(idx_x, idx_y) = STRONG_PIXEL;
				}
			}
		}

		c++; 
		if ( decider_2 <= 0 )
		{
			decider_2 += 2 * c + 1;
		}
		else
		{
			r-- ;
			decider_2 += 2 * ( c - r ) + 1;
		}

	}
}

double hough_detector::check_circle (const IMG_TYPE & hough, const int x, const int y , const int radius, kdtree & edgetree )
{
	const int width = hough.width ();
	const int height = hough.height ();

	int r = radius;
	int c = 0;
	int decider_2 = 1 - r ;

    struct kd_node_t * best = nullptr;
    double dist = 0.0;
    double area = 0.0; // discrete area 

	while (r >= c)
	{
		for ( int i = -1; i <= 1 ; i += 2)
		{
			for ( int j = -1; j <= 1; j += 2) 
			{
				int idx_x = x + r * i;
				int idx_y = y + c * j;
				
				if ( idx_x >= 0 && idx_x < width && idx_y >= 0 && idx_y < height) 
				{	
                    struct kd_node_t node = {{(double)idx_x, (double)idx_y}, nullptr, nullptr};
                    edgetree.nearest( &node, best, dist);
                    area += dist; 
				}
				idx_x = x + c * i;
				idx_y = y + r * j;
				
				if ( idx_x >= 0 && idx_x < width && idx_y >= 0 && idx_y < height) 
				{	
                    struct kd_node_t node = {{(double)idx_x, (double)idx_y}, nullptr, nullptr};
                    edgetree.nearest( &node, best, dist);
                    area += dist;
				}
			}
		}

		c++; 
		if ( decider_2 <= 0 )
		{
			decider_2 += 2 * c + 1;
		}
		else
		{
			r-- ;
			decider_2 += 2 * ( c - r ) + 1;
		}

	}
	return area; 
}

