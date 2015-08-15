#include "hough_include.h"

UINT nCr(UINT n, UINT r)
{
	if ( r > n )
		return 0.0; 
	if ( r == 0 || r == n )
		return 1.0;

	r = ( r > n - r)  ? n - r : r;
	UINT c = 1.0;
	for ( UINT i = 0; i < r ; i++ )
	{
		c = c * (n - i) / (i + 1);
	}
	return c;
}

void clamp( int & val, const int start, const int end ) 
{
    if ( val < start ) 
        val = start;
    else if ( val > end ) 
        val = end;
}
