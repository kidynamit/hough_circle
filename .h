#ifndef __UTILS_H__
#define __UTILS_H__

struct kd_node_t
{
    double x[MAX_DIM];
    struct kd_node_t *left, *right;
};

template<typename T>
class point
{
public:
    const T x;
    const T y; 

    point ( T x_ = 0, T y_ = 0) : x (x_), y(y_) {}
    ~point () {}
    
};

void clamp( int & val , const int start, const int end );

#endif
