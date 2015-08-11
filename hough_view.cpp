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
