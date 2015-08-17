#ifndef __HOUGH_MORPHOLOGY_H__
#define __HOUGH_MORPHOLOGY_H_

#include "hough_include.h"
#include "hough_strelement.h"

// erosion and dilation for feature extraction
// works with prescribed structuring elements
template <long unsigned int T>
IMG_TYPE morphology_operator( IMG_TYPE & hysterized_image, std::bitset<T> order, const std::vector<strelement> str_elements)
{
    const int width = hysterized_image.width () ;
    const int height = hysterized_image.height () ;
    
    IMG_TYPE out_image (width, height ) ;
    
    cimg_forXY(out_image, x, y) out_image(x , y) = hysterized_image(x, y);
    // erode out the points from the hysterized points
    int current_op = 0;
    
check:
    if ( current_op < order.size () )
    {
        if ( order[current_op] )
        {
            goto dilate;
        }
        else 
        {
            goto erode;
        }
    } 
    else 
    {
        goto finish;    
    }

// erode 

erode:
    current_op++;
    cimg_forXY(hysterized_image , row , col ) 
    {
        for ( int j = 0 ; j < str_elements.size (); j++ ) 
            for ( int i = 0; i < MAX_FEATURE_SIZE + 1; i++ )
            {
                if ( i == MAX_FEATURE_SIZE )
                {
                    out_image(row, col) = NULL_PIXEL;
                }
                else 
                {
                    // first  quad 
                    int r = row + str_elements[j].coords[i].x;
                    int c = col + str_elements[j].coords[i].y;

                    if ( r < width && r >= 0 && c >=0 && c < height ) 
                    {
                        if ( hysterized_image(r, c) == NULL_PIXEL) 
                            break;
                    }
                    //else { break ; } // for a stricter erosion
                }
            }
    }
    goto check;

dilate:
    current_op ++;
    cimg_forXY(hysterized_image , row , col ) 
    {
        for ( int j = 0 ; j < str_elements.size (); j++ ) 
            for ( int i = 0; i <MAX_FEATURE_SIZE; i++ )
            {
                int r = row + str_elements[j].coords[i].x;
                int c = col + str_elements[j].coords[i].y;

                if ( r < width && r >= 0 && c >=0 && c < height ) 
                {
                    if ( hysterized_image(row, col) == STRONG_PIXEL )
                        out_image(r, c) = STRONG_PIXEL; 
                }
                //else { break ; } // for a stricter erosion
            }
    }
    goto check;

finish:
    return out_image;
}


#endif
