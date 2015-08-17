#ifndef __HOUGH_INCLUDE_H__
#define __HOUGH_INCLUDE_H__

#define IMAGEPATH   "sandbox/images/"
#define PI          3.14159265359
#define MIN_RADIUS  5

#define LOG(TYPE, REASON) std::clog << #TYPE << "\t"<< __FILE__ << ":" << __LINE__ << "\t" << (REASON) << std::endl;

#define MIN(A, B) ( (A) > (B) ? (B) : (A) )
#define MAX(A, B) ( (A) > (B) ? (A) : (B) )

#include "CImg.h"
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <bitset>
#include <vector>

#include <omp.h>

typedef unsigned char UCHAR;
typedef unsigned int UINT;
typedef unsigned short USHORT;
typedef double PIXEL_TYPE;
using namespace cimg_library;

typedef CImg<PIXEL_TYPE> IMG_TYPE;
typedef CImgList<PIXEL_TYPE> IMG_LIST_TYPE;

// n choose r
UINT nCr( UINT n, UINT r );

void clamp( int & val , const int start, const int end );

static const PIXEL_TYPE STRONG_PIXEL = (PIXEL_TYPE)255;
static const PIXEL_TYPE WEAK_PIXEL = (PIXEL_TYPE)64;
static const PIXEL_TYPE NULL_PIXEL = (PIXEL_TYPE)0;

#define MAX_DIM 3

struct kd_node_t
{
    double x[MAX_DIM];
    struct kd_node_t *left, *right;
};

#endif
