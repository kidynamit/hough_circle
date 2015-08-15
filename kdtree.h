#ifndef __KD_TREE_H__
#define __KD_TREE_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
 
#include "hough_include.h"
class kdtree
{

	inline void swap(struct kd_node_t *x, struct kd_node_t *y) 
	{
		double tmp[MAX_DIM];
		memcpy(tmp,  x->x, sizeof(tmp));
		memcpy(x->x, y->x, sizeof(tmp));
		memcpy(y->x, tmp,  sizeof(tmp));
	}

	const UINT _dim;
	struct kd_node_t * _root;
public:

	kdtree(UINT dim=MAX_DIM);

	double
	dist(struct kd_node_t *a, struct kd_node_t *b);
	/* see quickselect method */
	struct kd_node_t*
	find_median(struct kd_node_t *start, struct kd_node_t *end, int idx);

	struct kd_node_t*
	make_tree(struct kd_node_t *t, int len, int i);

	void nearest(struct kd_node_t *root, struct kd_node_t *nd, int i,
			struct kd_node_t **best, double *best_dist);
};

#endif