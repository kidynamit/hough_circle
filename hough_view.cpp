#include "hough_detector.h"

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
		{ 
            if ( _min_threshold + diff < _max_threshold)
            { _min_threshold += diff; redetect = true; }
        }
		if (_display.is_keyK() )
		{ 
            if ( _min_threshold - diff > 0)
            { _min_threshold -= diff; redetect = true; }
        }
        if ( _display.is_keyU() )
		{ _max_threshold += diff; redetect = true; }
		if ( _display.is_keyI() )
        {
            if ( _max_threshold - diff > _min_threshold)
		    {  _max_threshold -= diff; redetect = true; }
        }
		

		if ( _display.is_keyV() )
		{ 
			{
				_hcd_threshold += 0.1;
				redetect_circle = true;
			}
		} 
		
		if ( _display.is_keyB() )
		{ 
			if ( _hcd_threshold - 0.1 > 0.0)
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

void hough_detector::print_usage() 
{
	std::clog << (
        "=============================================================================================\n"
		"\n                             ***HOUGH CIRCLE DETECTOR***\n"
		"                                                                               by WNJPAU001\n"
		            "CONTROLS:\n"
		"GAUSSIAN FILTER (optimal : ({ sigma : 2}, {window size : 5}))\n"
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        "       \"G\"               increase sigma value by 0.1\n"
        "       \"H\"               decrease sigma value by 0.1\n"
        "       \"T\"               increase window size by 2 ( limit 21 )\n"
        "       \"Y\"               decrease window size by 2 ( limit 3)\n"
        "\n"
        "CANNY EDGE DETECTOR (optimal: ({min. threshold : 100.0 }, {max. threshold : 200.0}\n"
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        "       \"J\"               increase min. hysteresis threshold by 10 (limit max. threshold)\n"
        "       \"K\"               decrease min. hysteresis threshold by 10 (limit 0)\n"
        "       \"U\"               increase max. hysteresis threshold by 10\n"
        "       \"I\"               decrease max. hysteresis threshold by 10 (limit min. threshold)\n"
        "\n"
        "HOUGH CIRCLE DETECTOR (optimal :( { high threshold : 1.1} )\n"
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        "       \"V\"               increase high threshold by 0.1 \n"
        "       \"B\"               decrease high threshold by 0.1 (limit 0.0) \n"
        "\n NOTE: Title bar acts as status bar. Controls only work in view of the display.\n"
        "=============================================================================================\n"
		);
}

