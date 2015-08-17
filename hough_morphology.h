#ifndef __HOUGH_MORPHOLOGY_H__
#define __HOUGH_MORPHOLOGY_H_

#include "hough_include.h"

// erosion for feature extraction
template <long unsigned int T>
IMG_TYPE morphology_operator( IMG_TYPE & hysterized_image, std::bitset<T> order)
{
    const int width = hysterized_image.width () ;
    const int height = hysterized_image.height () ;
    
    IMG_TYPE out_image (width, height ) ;
    
    const int strel_num = 3; // number of structuring elements
    // defining structuring elements for a small circle of radius 1
    const int num_strels = 12;
    
    int struct_elements [num_strels][strel_num][2] = 
                    { 
                        { {-1, 0}, {0, 0}, {1, 0} }, // STRAIGHT LINES
                        { {0, -1}, {0, 0}, {0, 1} }, 
                        { {1, 1}, {0, 0}, {1, 1} }, 
                        { {1, -1}, {0, 0}, {-1, 1} },     
                        { {-1, 0}, {0, 0}, {1, 1} }, // ARCS
                        { {-1, 0}, {0, 0}, {1, -1} }, 
                        { {1, 0}, {0, 0}, {-1, 1} }, 
                        { {1, 0}, {0, 0}, {-1, -1} }, 
                        { {0, -1}, {0, 0}, {-1, 1} }, 
                        { {0, -1}, {0, 0}, {1, 1} }, 
                        { {0, 1}, {0, 0}, {1, -1} }, 
                        { {0, 1}, {0, 0}, {-1, -1} }
                        /*   // corners
                        { {0, 0}, { -1, 1}, {-1, -1} }, 
                        { {0, 0}, { -1, 1}, {1, 1} },
                        { {0, 0}, {1 ,1} , {1, -1} }, 
                        { {0, 0}, {1, -1}, {-1, -1}}
                        */
                    };
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
        for ( int j = 0 ; j < num_strels; j++ ) 
            for ( int i = 0; i < strel_num + 1; i++ )
            {
                if ( i == strel_num )
                {
                    out_image(row, col) = NULL_PIXEL;
                }
                else 
                {
                    // first  quad 
                    int r = row + struct_elements[j][i][0];
                    int c = col + struct_elements[j][i][1];

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
        for ( int j = 0 ; j < num_strels; j++ ) 
            for ( int i = 0; i < strel_num + 1; i++ )
            {
                int r = row + struct_elements[j][i][0];
                int c = col + struct_elements[j][i][1];

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
