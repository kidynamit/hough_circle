#include "utils.h"

void clamp( int & val, const int start, const int end ) 
{
    if ( val < start ) 
        val = start;
    else if ( val > end ) 
        val = end;
}
