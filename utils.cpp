#include "utils.h"

/**
 * clamps a value between a start and end value
 */
void clamp( int & val, const int start, const int end ) 
{
    if ( val < start ) 
        val = start;
    else if ( val > end ) 
        val = end;
}
