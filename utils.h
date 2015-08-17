#ifndef __UTILS_H__
#define __UTILS_H__

#define MAX_DIM 3

struct kd_node_t
{
    double x[MAX_DIM];
    struct kd_node_t *left, *right;
};


struct point
{
	const int x;
    const int y;   
};

void clamp( int & val , const int start, const int end );

#endif
