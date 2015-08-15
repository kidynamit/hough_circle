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

		for ( int i = 0; i < 4 ; i++)
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
    
    double GX[width][height], GY[width][height], THETA[width][height];

	const int SobelY [3][3] =
		{ 
			{ -1, -2, -1}, 
			{0, 0, 0}, 
			{1, 2, 1}
		};
    const int Vertical [3][3] =
        {
            { -1, 2, -1}, 
            { -1, 2, -1}, 
            { -1, 2, -1}
        };
    
    const int Horizontal [3][3] =
        {
            { -1, -1, -1}, 
            { 2, 2, 2}, 
            { -1, -1, -1}
        };

    const int LeftDiagonal [3][3] =
        {
            { 2, -1, -1}, 
            { -1, 2, -1}, 
            { -1, -1, 2}
        };

    const int RightDiagonal [3][3] =
        {
            { -1, -1, 2}, 
            { -1, 2, -1}, 
            { 2, -1, -1}
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
			// calculate anglea
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
            if (theta < 0.0 )
                theta = theta + PI;

            THETA[row][col] = PI/2.0 - (theta) + 2.0 * PI; 

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
				int i = -1;
				// loop through the neighbours and activate the weak pixels
				for ( ; i <= 1 ; i++ )
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
                            break;
						}
					}
				}
			    if ( i > 1 )
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
   

    cimg_forXY(hysterized_image, row, col)
    {
        if (hysterized_image(row, col) != STRONG_PIXEL )
            output_image(row, col) = NULL_PIXEL;
        else
        {
            
            double vert_sum = 0 , hor_sum = 0, leftdiag_sum = 0, rightdiag_sum = 0;
            for ( int i = -1; i<= 1; i++ )
            {
                for ( int j = -1; j <= 1; j++ )
                {
                    int r = row + i;
                    int c = col + j;
                    
                    if ( r < width && r >= 0 && c < height && c >= 0)
                    {
                        vert_sum += Vertical[i + 1][j + 1] * THETA[r][c] / THETA[row][col];    
                        hor_sum += Horizontal[i + 1][j + 1] * THETA[r][c]/ THETA[row][col];    
                        leftdiag_sum += LeftDiagonal[i + 1][j + 1] * THETA[r][c]/ THETA[row][col];    
                        rightdiag_sum += RightDiagonal[i + 1][j + 1] * THETA[r][c]/ THETA[row][col];    
                    }

                }
            }
            double max_sum = vert_sum ;
            if (max_sum < hor_sum ) max_sum = hor_sum;
            if (max_sum < leftdiag_sum ) max_sum = leftdiag_sum;
            if (max_sum < rightdiag_sum ) max_sum = rightdiag_sum;
        
            double thresh = 0.02 *THETA[row][col];
            if ( max_sum > thresh)
            {
                output_image(row, col) = NULL_PIXEL;
            }
            else
            {
                output_image(row, col) = STRONG_PIXEL;
            }
        }
    }

	//std::clog << "max " << max_magnitude << " min "<< min_magnitude << std::endl;
	_images[2] = output_image;  
    
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
        int num_pixels = 0;
        int total = 0;
		for ( UINT x = 0 ;  x < width ; x ++ ) 
        {
            for ( UINT y = 0; y < height ; y++ )
		    {
			    if ( _images[2](x, y) == STRONG_PIXEL)
			    {            
                    num_pixels += accumulate_circle ( hough, x, y, radius);
                    total ++;
			    }
		    }
        }

        double threshold =  (num_pixels )/ (total * _hcd_threshold); 
        // double threshold = _hcd_threshold ;
		for ( UINT x = 0;  x < width ; x ++ ) 
        {
            for ( UINT y = 0; y < height ; y++ )
		    {
			    if ( ( hough(x, y)) > threshold)
			    {
				    draw_circle ( out_image, x, y, radius);
			    } 
		    }
        }
	}
    _images[3] = hough;
    _images[4] = out_image;
    
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
