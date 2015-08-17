#include "hough_detector.h"


hough_detector::hough_detector (const char * file, double sigma , int gaussian_window, double max_thresh, double min_thresh)
	: _sigma(sigma), 
	_gaussian_window(gaussian_window), 
	_max_threshold(max_thresh), 
	_min_threshold(min_thresh)
{
	_filename = const_cast<char*>(std::string(file).c_str());
	_hcd_threshold = 1.1;   
	_hcd_low_threshold = 0.85 ;
	init ();
}

hough_detector::~hough_detector()
{

}


void hough_detector::init()
{
	print_usage ();
	_display = CImgDisplay( _images, DISPLAY_TITLE, 0);
    _display.set_title ( (DISPLAY_TITLE " ( Please view instructions on the console ) ") );
    cimg::wait (500) ;
    _display.set_title ( DISPLAY_TITLE);
	std::string s (_filename);
	try
	{
		_images.push_back( IMG_TYPE(_filename).normalize(0, 255) );
		
		const int width = _images[0].width();
		const int height = _images[0].height();
		
		_min_radius = MIN_RADIUS;
		_max_radius = ( width > height ) ? height / 2 : width / 2;

		for ( int i = 0; i < 4; i++)
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
			for ( int i = -start ; i <= start ; i++ )
			{
				for ( int j = -start ; j <= start ; j++ )
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
    IMG_TYPE canny_image (width, height );	
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
				x = 1 ; y = 1;
			}
			else if ( theta >= -3.0 * PI/8.0 && theta < -PI/8.0 )
			{
				x = 1; y = -1;
			}
			else if ( theta > -3.0 *PI/8.0 || theta < -3.0 * PI/8.0 )
			{
				x = 1; y = 0;
			}
			else 
			{ LOG(INFO, "bug");}

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
		    canny_image (row, col) = magnitude;	
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

	_edge_pixels.clear();
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
	

    // perform operation 
    // set the morphology operations order
    const long unsigned int depth = 40 ;
    std::bitset<depth> closing, opening; 
    for ( size_t i = 0 ; i < depth ; i++ )
    {
        closing [i] = (i + 1) % 2;
        opening [i] = i % 2;
    }

    //for ( int i = depth/2; i < depth; i++ )
        //order[i] = (i + 1) % 2
    _images[3] = morphology_operator (hysterized_image, opening, Corners);
    _images[3] = morphology_operator (_images[3], closing, Arcs);
    // record the edge pixels
    cimg_forXY(_images[3], row, col ) 
    {
        if (_images[3](row, col) == STRONG_PIXEL)
		{
		    _edge_pixels.push_back({{(double)row, (double)col},  nullptr, nullptr});
		}
    }


    _display.set_title ( (DISPLAY_TITLE " ... Ready!") );
	cimg::wait (1000);
	RERENDER;
	_display.set_title ( (DISPLAY_TITLE) );
}


