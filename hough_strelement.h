#ifndef __HOUGH_STRUCTURING_ELEMENT_H__
#define __HOUGH_STRUCTURING_ELEMENT_H__

#define MAX_FEATURE_SIZE 3

#include "hough_include.h"

struct strelement
{
  	const point coords [MAX_FEATURE_SIZE];
};

const static strelement __arcs[8] =
{
    { { {-1, 0}, {0, 0}, {1, -1} }  } , 
    { { {-1, 0}, {0, 0}, {1, 1}  }  } , // ARCS
    { { {1, 0}, {0, 0}, {-1, 1} }  } , 
    { { {1, 0}, {0, 0}, {-1, -1} }  } , 
    { { {0, -1}, {0, 0}, {-1, 1} }  } , 
    { { {0, -1}, {0, 0}, {1, 1} }  } , 
    { { {0, 1}, {0, 0}, {1, -1} }  } , 
    { { {0, 1}, {0, 0}, {-1, -1} }  } 
};

const static strelement __straight_lines[4] =
{
	{ { {-1, 0}, {0, 0}, {1, 0} }  } , // STRAIGHT LINES
	{ { {0, -1}, {0, 0}, {0, 1} }  } , 
	{ { {1, 1}, {0, 0}, {1, 1} }  } , 
	{ { {1, -1}, {0, 0}, {-1, 1} }  } 
};

const static strelement __corners[4] = 
{
	{ { {0, 0}, {-1, 1}, {-1, -1} }  } , 
	{ { {0, 0}, {-1, 1}, {1, 1} }  } ,
	{ { {0, 0}, {1, 1} , {1, -1} }  } , 
	{ { {0, 0}, {1, -1}, {-1, -1} }  } 
};

const static std::vector<strelement>  Arcs ( __arcs , __arcs + 8);
const static std::vector<strelement>  StraightLines ( __straight_lines, __straight_lines + 4);
const static std::vector<strelement>  Corners ( __corners, __corners + 4);
#endif
