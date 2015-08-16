#include "hough_detector.h"

#define RERENDER { _display.resize ( _display, false).display(_images); }//_display.show ();}

hough_detector::hough_detector (const char * file, double sigma , int gaussian_window, double max_thresh, double min_thresh)
	: _sigma(sigma), 
	_gaussian_window(gaussian_window), 
	_max_threshold(max_thresh), 
	_min_threshold(min_thresh)
{
	_filename = const_cast<char*>(std::string(file).c_str());
	_hcd_threshold = 2.4;   
	_hld_threshold = 2.4;
	init ();
}

hough_detector::~hough_detector()
{

}

void hough_detector::print_usage() 
{
	std::clog << (
		"***Hough Circle Detector***\n"
		"                           by kidynamit\n"
		"CONTROLS:\n"
		"   G + mouse-wheel     alter the pre-filtering gaussian sigma value\n"
		);
}

void hough_detector::init()
{
	print_usage ();
	std::string s (_filename);
	try
	{
		_images.push_back( IMG_TYPE(_filename).normalize(0, 255) );
		
		const int width = _images[0].width();
		const int height = _images[0].height();
		
		_min_radius = MIN_RADIUS;
		_max_radius = ( width > height ) ? height / 2 : width / 2;

		for ( int i = 0; i < 5; i++)
			_images.push_back( IMG_TYPE(width, height));

		std::stringstream ss ;
		int spectrum = _images[0].spectrum ();
		ss <<  "image \"" << s << "\" successfully loaded \t[" << width << " x " << height << "] " << spectrum ;  
		LOG(INFO, ss.str().c_str());
	}
	catch (...)
	{
		LOG(ERROR, "unable to read file \"" + s + "\"");
	}
	double ** kernel = nullptr;
	double kernel_sum = 0.0;
	
	_display = CImgDisplay( _images, DISPLAY_TITLE, 0);
	gaussian_filter ( kernel, kernel_sum ); 
	canny_edge_detector () ;
	hough_circle_detector ();
	if ( kernel )
	{
		for (int i = 0; i < _gaussian_window ; i ++)
			delete [] kernel [i]; 
		delete [] kernel;
	}
	
}

void hough_detector::event_loop()
{
	double ** kernel = nullptr;
	double kernel_sum = 0.0;
	while (!_display.is_closed () && !_display.is_keyESC())
	{
		bool reblur = false;
		bool redetect = false;
		bool redetect_circle = false;
		if ( _display.is_keyG() )//&& _display.wheel()) 
		{ 
			_sigma += 0.1 ;
			//_gaussian_window += 2 * ( h ? -1 : 1);
			reblur = true;
		}       
		
		if ( _display.is_keyH()) //&& _display.wheel()) 
		{
			
			//_gaussian_window += 2 * _display.wheel();
			 //_gaussian_window += 2 * ( h ? -1 : 1);
			_sigma -= 0.1; 
			reblur = true;
		}
		if ( _display.is_keyT()) 
		{
			_gaussian_window += 2;
			if ( _gaussian_window / 2 > 10) 
				_gaussian_window = 21;
			else 
				reblur = true;
		}
		if (_display.is_keyY() )
		{
			_gaussian_window -= 2;
			if ( _gaussian_window / 2 < 1 )
				_gaussian_window = 3;
			else 
				reblur = true; 
		}
		int diff = 10;
		if (_display.is_keyJ() )
		{ _min_threshold += diff; redetect = true; }
		if (_display.is_keyK() )
		{ _min_threshold -= diff; redetect = true; }
		if ( _display.is_keyU() )
		{ _max_threshold += diff; redetect = true; }
		if ( _display.is_keyI() )
		{ _max_threshold -= diff; redetect = true; }

		if ( _display.is_keyO() )
		{ _min_radius += 2; redetect = true; }
		if ( _display.is_keyP() )
		{ 
			if ( _min_radius / 2 != 0)
			{
				_min_radius -= 2;
				redetect = true;
			}       
		} 
		if ( _display.is_keyV() )
		{ 
			if ( _hcd_threshold + 0.1 < STRONG_PIXEL )
			{
				_hcd_threshold += 0.1;
				redetect_circle = true;
			}
		} 
		
		if ( _display.is_keyB() )
		{ 
			if ( _hcd_threshold - 0.1 != 0.0)
			{
				_hcd_threshold-= 0.1;
				redetect_circle = true;
			}
		} 

		if ( _display.is_keyQ() )
		{ 
			if ( _hld_threshold + 0.1 < STRONG_PIXEL )
			{
				_hld_threshold += 0.1;
				redetect = true;
			}       
		} 
		
		if ( _display.is_keyW() )
		{ 
			if ( _hld_threshold - 0.1 != 0.0)
			{
				_hld_threshold -= 0.1;
				redetect = true;
			}
		} 
		if ( reblur )
		{ 
			kernel = nullptr;
			kernel_sum = 0.0;
			gaussian_filter (kernel, kernel_sum);
			redetect = true;
			redetect_circle = true;
		}
		if (redetect)
		{
			canny_edge_detector () ;
			redetect_circle = true;
		}

		if (redetect_circle)
		{
			hough_circle_detector () ;
		}
	}
}

void hough_detector::gaussian_filter( double **& kernel, double& kernel_sum, const int idx_in, const int idx_out)
{
	std::stringstream ss;
	ss << "Applying Gaussian filter: \t\t\t" << _sigma << " , "<< _gaussian_window; 
	LOG(INFO, ss.str().c_str() );
	_display.set_title ( (DISPLAY_TITLE " ... Applying Gaussian Filter") );
	
	if (_gaussian_window %2 == 0)
	{
		LOG(ERROR, "gaussian window size is even");
		return;
	}
	int width = _images[idx_in].width () , height = _images[idx_in].height(); 
	IMG_TYPE output_image (width, height );

	int start = _gaussian_window / 2;
	double variance = _sigma * _sigma; 
	
	
	
	{
		kernel = new double * [_gaussian_window];
		for ( int i = 0 ; i < _gaussian_window ; i++ )
			kernel[i] = new double[_gaussian_window];
		kernel_sum = 0.0;
		for (int i = -start; i <=start; i++ )
		{
			for ( int j = -start; j <= start ; j++ )
			{
				double r = double(j) * double(j) + double(i) * double(i);
				kernel[i + start][j + start ] = (1 / ( 2.0 * PI * variance )) * (exp ( (-r) / (2.0 * variance)) );
				//kernel[i ][j ] = (nCr(_gaussian_window - 1, i) * nCr(_gaussian_window - 1, j));
				kernel_sum += kernel[i + start][j + start];		
			}
		}
	}   
	
	for ( int row = 0 ; row < width ; row++ )
	{
		for ( int col = 0 ; col < height ; col++ )
		{
			double sum = 0.0;
			// weighted sum of the kernel
			for ( int i = start ; i >= -start ; i-- )
			{
				for ( int j = start ; j >= -start ; j-- )
				{
					int r = row + i;
					int c = col + j;

					clamp (r, 0, width - 1);
					clamp (c, 0, height - 1);

					sum = sum + (kernel[i + start][j + start] * _images[idx_in](r, c)) ;/// kernel_sum;
				}
			}
			output_image(row, col) = (PIXEL_TYPE)sum;
		}
	}
	
	_images[idx_out] = (output_image);

	for ( int i = 0 ; i < _gaussian_window ; i++ )
		delete [] kernel[i];
	delete [] kernel ;
	kernel = nullptr;
	
	_display.set_title ( (DISPLAY_TITLE " ... Ready!") );
	cimg::wait (1000);
	RERENDER;
	_display.set_title ( (DISPLAY_TITLE) );
}

void hough_detector::canny_edge_detector () 
{
	std::stringstream ss;
	ss << "Applying Canny Edge Detection: \t\t\t"  <<  _min_threshold << " , "<< _max_threshold; 
	LOG(INFO, ss.str().c_str() );
	
	//_display.close();
	const int width = _images[0].width ();
	const int height = _images[0].height ();

	// apply sobel filter
	IMG_TYPE output_image ( width , height) ;
	IMG_TYPE hysterized_image (width, height);
	
	double GX[width][height], GY[width][height];

	const int SobelY [3][3] =
		{ 
			{ -1, -2, -1}, 
			{0, 0, 0}, 
			{1, 2, 1}
		};
	// convolve image 
	for ( int row = 0 ; row < width ; row++ )
	{
		for ( int col = 0 ; col < height ; col++ )
		{
			double sumX = 0.0, sumY = 0.0;
			for ( int i = -1; i <= 1 ; i ++ )
			{
				for ( int j = -1 ; j <= 1 ; j++ )
				{
					int r = row + j;
					int c = col + i;

					clamp (r, 0, width - 1);
					clamp (c, 0, height - 1);

					{
						sumX += _images[1](r, c) * SobelY[j + 1][i + 1];
						sumY += _images[1](r, c) * SobelY[i + 1][j + 1];
					}
				}
			}
			GX[row][col] = sumX;
			GY[row][col] = sumY;
		}    
	}

	double max_magnitude = -1.0, min_magnitude = -1.0;
	
	// non-maximum suppression
	for ( int row = 0 ; row < width ; row++ )
	{
		for ( int col = 0 ; col < height ; col++ ) 
		{
			int x = 0, y = 0; // offset indicies of the position to check for
			// calculate angle
			double theta;
			if ( GX[row][col] != 0.0 ) 
				theta = atan(GY[row][col] / GX[row][col]); 
			else
				theta = PI/2.0;
			// classify
			if ( theta >= -PI/8.0 && theta < PI/8.0 ) 
			{
				x = 0; y = 1;
			}
			else if ( theta >= PI/8.0 && theta < 3.0 * PI/8.0 )
			{ 
				x = 1 ; y = -1;
			}
			else if ( theta >= -3.0 * PI/8.0 && theta < -PI/8.0 )
			{
				x = 1; y = 1;
			}
			else if ( theta > -3.0 *PI/8.0 || theta < -3.0 * PI/8.0 )
			{
				x = 1; y = 0;
			}
			else 
			{}

			int r_left = row - y; 
			int r_right = row + y;
			
			int c_left = col - x;
			int c_right = col + x ;

			clamp ( r_left, 0, width - 1);
			clamp ( r_right, 0, width - 1);
			
			clamp ( c_left, 0, height - 1);
			clamp ( c_right, 0, height - 1);

			double magnitude = (double)hypot(GX[row][col], GY[row][col]);
			if ( min_magnitude == -1.0 || min_magnitude > magnitude)
				min_magnitude = magnitude;

			if ( magnitude > max_magnitude || max_magnitude == -1.0)
				max_magnitude = magnitude;

			double magnitude_left = (double)hypot(GX[r_left][c_left], GY[r_left][ c_left]);
			double magnitude_right = (double)hypot(GX[r_right][c_right], GY[r_right][ c_right]);
			
			if ( magnitude > magnitude_left && magnitude > magnitude_right && magnitude > _min_threshold)
			{
				if ( magnitude > _max_threshold )
					output_image(row, col) = STRONG_PIXEL;
				else 
					output_image(row, col) = WEAK_PIXEL;
			}
			else
				output_image(row, col) = NULL_PIXEL;
			hysterized_image(row, col) = NULL_PIXEL;
		}
				
	}
	
	for ( int row = 0; row < width; row++ )
	{
		for ( int col = 0 ; col < height ; col++ )
		{
			if ( output_image(row, col) == WEAK_PIXEL)
			{
				bool found = false;
				// loop through the neighbours and activate the weak pixels
				for ( int i = -1 ; i <= 1 && !found; i++ )
				{
					for ( int j = -1 ; j <= 1 ; j++ )
					{
						int r = row + i;
						int c = col + j;
						
						clamp (r, 0, width - 1);
						clamp (c, 0, height - 1);
						if ( output_image(r, c) == STRONG_PIXEL) 
						{
							hysterized_image(row, col) = STRONG_PIXEL;
							found = true;
							break;
						}
					}
				}
				if ( !found )
				{
					hysterized_image(row, col) = NULL_PIXEL;
				}
			}
			else 
			{
				hysterized_image(row, col) = output_image(row, col);
			}
			
		}
	}
	//hough_line_detector(hysterized_image);   
	//std::clog << "max " << max_magnitude << " min "<< min_magnitude << std::endl;
	_images[2] = hysterized_image;

	
	_display.set_title ( (DISPLAY_TITLE " ... Ready!") );
	cimg::wait (1000);
	RERENDER;
	_display.set_title ( (DISPLAY_TITLE) );
}

void hough_detector::hough_circle_detector ()
{
	std::stringstream ss;
	ss << "Applying Hough Circle Detection: \t\t" << _hcd_threshold <<" - "  <<  _min_radius<< " , "<< _max_radius; 
	LOG(INFO, ss.str().c_str() );

	_display.set_title ( (DISPLAY_TITLE " ... Applying Hough Circle Detection") );
	
	// _display.close();
	const int width = _images[0].width();
	const int height = _images[0].height();
	
	// construct list 
	IMG_TYPE out_image (width, height);
	IMG_TYPE accumulator (width, height );
	cimg_forXY(out_image, x, y) out_image(x, y) = NULL_PIXEL;

	IMG_TYPE hough (width, height);
	for ( UINT radius = _min_radius; radius < _max_radius ; radius ++ ) 
	{
		
		cimg_forXY(hough, x, y) hough(x, y) = NULL_PIXEL;
		int cum_circumference = 0;
		int total = 0;
		for ( int x = 0 ;  x < width ; x ++ ) 
		{
			for ( int y = 0; y < height ; y++ )
			{
				if ( _images[2](x, y) == STRONG_PIXEL)
				{            
					cum_circumference += accumulate_circle ( hough, x, y, radius);
					total ++;
				}
			}
		}

		double threshold =  (cum_circumference )/ (total * _hcd_threshold); 
		// double threshold = _hcd_threshold ;
		for ( int x = 0;  x < width ; x ++ ) 
		{
			for ( int y = 0; y < height ; y++ )
			{
				if ( ( hough(x, y)) > threshold)
				{
					draw_circle ( out_image, x, y, radius);
				} 
			}
		}
	}
	_images[4] = hough;
	_images[5] = out_image;
	
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

void hough_detector::draw_circle (IMG_TYPE & hough, const int x, const int y , const int radius )
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

void hough_detector::hough_line_detector (IMG_TYPE & in_image)
{
	std::stringstream ss;
	ss << "Applying Hough Line Detection: \t\t" << _hld_threshold;
	LOG(INFO, ss.str().c_str() );

	_display.set_title ( (DISPLAY_TITLE " ... Applying Hough Line Detection") );
	
	// _display.close();
	const int width = _images[0].width();
	const int height = _images[0].height();
	
	// construct list

	IMG_TYPE hough (width, height);

	for ( int x = 0 ;  x < width ; x ++ ) 
	{
		for ( int y = 0; y < height ; y++ )
		{
			if ( in_image(x, y) == STRONG_PIXEL)
			{            
				accumulate_line ( hough, x, y);
			}
		}
	}

	double threshold =  _hld_threshold; 
	// double threshold = _hcd_threshold ;
	for ( int x = 0;  x < width ; x ++ ) 
	{
		for ( int y = 0; y < height ; y++ )
		{
			if ( ( hough(x, y)) > threshold)
			{
				unmark_line ( in_image, x, y);
			} 
		}
	}
	_images[3] = hough;
	
	_display.set_title ( (DISPLAY_TITLE " ... Ready!") );
	cimg::wait (1000);
	RERENDER;
	_display.set_title ( (DISPLAY_TITLE) );
}

void hough_detector::accumulate_line( IMG_TYPE & hough, const int m, const int c)
{
	const int width = hough.width (); const int height = hough.height ();
	int x1 = 0; 
	int y1 = c;
	clamp (y1, 0, height);

	int x2 = width ;
	int y2 = m * x2 + c;
	clamp (y2, 0, height);

	// Bresenham's line algorithm
	const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
	if(steep)
	{
		std::swap(x1, y1);
		std::swap(x2, y2);
	}

	if(x1 > x2)
	{
		std::swap(x1, x2);
		std::swap(y1, y2);
	}

	const int dx = x2 - x1;
	const int dy = fabs(y2 - y1);

	float error = dx / 2.0f;
	const int ystep = (y1 < y2) ? 1 : -1;
	int y = (int)y1;

	const int maxX = (int)x2;

	for(int x=(int)x1; x<maxX; x++)
	{
		if(steep)
		{
			if ( y < width && y >= 0 && x < height && x >= 0)
			{
				hough(y, x) ++;
			} 
		}
		else
		{
			if ( x < width && x >= 0 && y < height && y >= 0)
			{
				hough(x, y) ++;
			} 
		}
 
		error -= dy;
		if(error < 0)
		{
			y += ystep;
			error += dx;
		}
	}
}

void hough_detector::unmark_line( IMG_TYPE & hough, const int m, const int c)
{
	const int width = hough.width (); const int height = hough.height ();

	int x1 = 0; 
	int y1 = c;
	clamp (y1, 0, height);

	int x2 = width ;
	int y2 = m * x2 + c;
	clamp (y2, 0, height);

	// Bresenham's line algorithm
	const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
	if(steep)
	{
		std::swap(x1, y1);
		std::swap(x2, y2);
	}

	if(x1 > x2)
	{
		std::swap(x1, x2);
		std::swap(y1, y2);
	}

	const int dx = x2 - x1;
	const int dy = fabs(y2 - y1);

	float error = dx / 2.0f;
	const int ystep = (y1 < y2) ? 1 : -1;
	int y = (int)y1;

	const int maxX = (int)x2;

	for(int x=(int)x1; x<maxX; x++)
	{
		if(steep)
		{
			if ( y < width && y >= 0 && x < height && x >= 0)
			{
				hough(y, x) = NULL_PIXEL;
			} 
		}
		else
		{
			if ( x < width && x >= 0 && y < height && y >= 0)
			{
				hough(x, y) = NULL_PIXEL;
			} 
		}
 
		error -= dy;
		if(error < 0)
		{
			y += ystep;
			error += dx;
		}
	}
}

