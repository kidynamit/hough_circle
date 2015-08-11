#include "hough_detector.h"

hough_detector::hough_detector (const char * file, double sigma , UINT gaussian_window)
	: _sigma(sigma), 
	_gaussian_window(gaussian_window)
{
	_filename = const_cast<char*>(std::string(file).c_str());
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
		_image = IMG_TYPE(_filename);
		std::stringstream ss ;
		ss <<  "image \"" << s << "\" successfully loaded \t[" << _image.width() << " x " << _image.height() << "] " << _image.spectrum() ;  
		LOG(INFORM, ss.str().c_str());
	}
	catch (...)
	{
		LOG(ERROR, "unable to read file \"" + s + "\"");
	}
	_display = CImgDisplay( _image, "Hough circle detector", 0);
	gaussian_filter ();
	//_image.blur(_sigma, false, true);
}

void hough_detector::event_loop()
{
	while (!_display.is_closed () && !_display.is_keyESC())
	{
		if ( _display.wheel() && _display.is_keyG() ) 
		{
			_sigma += 2 * _display.wheel (); 
			// _image.blur(_sigma, false, true);
			gaussian_filter ();
		}
	}
}

void hough_detector::gaussian_filter()
{
	if (_gaussian_window %2 == 0)
	{
		LOG(ERROR, "gaussian window size is even");
		return;
	}
	
	IMG_TYPE output_image (_image.width() ,_image.height());

	int start = _gaussian_window / 2;
	double variance = _sigma * _sigma; 
	// construct the gaussian kernel
	double ** kernel = new double * [_gaussian_window];
	for ( UINT i = 0 ; i < _gaussian_window ; i++ )
		kernel[i] = new double[_gaussian_window];
	double sum = 0.0;
	for ( UINT i = 0; i < _gaussian_window; i++ )
	{
		for ( UINT j = 0; j < _gaussian_window ; j++ )
		{
			//double r = double(j) * double(j) + double(i) * double(i);
			//kernel[i + start][j + start ] = (1 / ( 2.0 * PI * variance )) * (exp ( (-r) / (2.0 * variance)) );
			kernel[i ][j ] = (nCr(_gaussian_window - 1, i) * nCr(_gaussian_window - 1, j));
			sum += kernel[i ][j ];		
		}
	}
	for ( int i = -start; i <= start; i++ )
	{
		for ( int j = -start; j <= start ; j++ )
		{
			kernel[i + start][j + start] = kernel[i + start][j + start] / sum;
		}
	}



	for ( int row = 0 ; row < _image.width() ; row++ )
	{
		for ( int col = 0 ; col < _image.height () ; col++ )
		{
			double sum = 0.0;
			// weighted sum of the kernel
			for ( int i = start ; i >= -start ; i-- )
			{
				for ( int j = start ; j >= -start ; j-- )
				{
					UINT r, c;
					if ( row + i < 0 )
						r = 0;
					else if ( row + i > _image.height() - 1)
						r = _image.height() - 1;
					else 
						r = (row + i) ;

					if ( col + i < 0 )
						c = 0;
					else if ( row + i > _image.width() - 1 ) 
						c = _image.height() - 1;
					else
						c = (col + i); 

					sum = sum + (kernel[i + start][j + start] * _image(r, c));
				}
			}
			output_image(row, col) = (PIXEL_TYPE)sum;// / (float(_gaussian_window) * float(_gaussian_window) ));
			//output_image(row, col ) = _image(row, col);
			//std::clog << _image(row, col) << " ";
		}
		//std::clog << std::endl;
	}

	for ( UINT i = 0 ; i < _gaussian_window ; i++ )
		delete [] kernel[i];
	delete [] kernel ;
	// _gaussian_display = CImgDisplay (output_image, "gaussian filter");
}

