#ifndef __HOUGH_INCLUDE_H__
#define __HOUGH_INCLUDE_H__

#define IMAGEPATH   "sandbox/images/"
#define PI          3.14159265359

#define LOG(TYPE, REASON) std::clog << #TYPE << "\t"<< __FILE__ << ":" << __LINE__ << "\t" << (REASON) << std::endl;

#include "CImg.h"
#include <iostream>
#include <string>
#include <sstream>

typedef unsigned char UCHAR;
typedef unsigned int UINT;
typedef unsigned short USHORT;
typedef UCHAR PIXEL_TYPE;
using namespace cimg_library;

typedef CImg<PIXEL_TYPE> IMG_TYPE;

// n choose r
UINT nCr( UINT n, UINT r );

#endif
