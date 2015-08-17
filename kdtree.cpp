#include "kdtree.h"

/**
 * recursively finds the nearest node and the distance to the nearest 
 * node in the kd tree
 */
void kdtree::nearest(struct kd_node_t *& root, struct kd_node_t *& nd, int i,
		struct kd_node_t *&best, double &best_dist)
{
	double d, dx, dx2;
 
	if (!root) return;
	d = dist(root, nd);
	dx = root->x[i] - nd->x[i];
	dx2 = dx * dx;

	if (!best || d < best_dist) 
	{
		best_dist = d;
		best = root;
	}
 
	/* if chance of exact match is high */
	if (!best_dist) return;
 
	if (++i >= (int)_dim) i = 0;
 
	nearest(dx > 0 ? root->left : root->right, nd, i, best, best_dist);
	if (dx2 >= best_dist) return;
	nearest(dx > 0 ? root->right : root->left, nd, i, best, best_dist);
}

/**
 * constructs a kd tree given an array of nodes
 */
struct kd_node_t*
kdtree::make_tree(struct kd_node_t * t, int len, int i)
{
	struct kd_node_t *n;
 
	if (!len) return 0;
 
	if ((n = find_median(t, t + len, i))) 
    {
		i = (i + 1) % _dim;
		n->left  = make_tree(t, n - t, i);
		n->right = make_tree(n + 1, t + len - (n + 1), i);
	}
	return n;
}

/**
 * calculates the median between two nodes
 * uses the iterative quick select method
 */
struct kd_node_t*
kdtree::find_median(struct kd_node_t *& start, struct kd_node_t * end, int idx)
{
	if (end <= start) return NULL;
	if (end == start + 1)
		return start;
 
	struct kd_node_t *p, *store, *md = start + (end - start) / 2;
	double pivot;
	while (1) 
    {
		pivot = md->x[idx];
        struct kd_node_t * prev = end - 1;
		swap_node(md, prev);
		for (store = p = start; p < end; p++) {
			if (p->x[idx] < pivot) {
				if (p != store)
					swap_node(p, store);
				store++;
			}
		}
		swap_node(store, end - 1);
 
		/* median has duplicate values */
		if (store->x[idx] == md->x[idx])
			return md;
 
		if (store > md) end = store;
		else        start = store;
	}
}
/**
 * wraps the nearest function to work with the defined root node
 */
void kdtree::nearest(struct kd_node_t * node, struct kd_node_t *& best, double & dist )
{
    if ( _root ) 
    {
        nearest (_root , node, 0, best, dist);
    }
}
/**
 * returns the square of the Euclidean distance of two points. 
 * avoids floating points
 */
double
kdtree::dist(struct kd_node_t *& a, struct kd_node_t *& b)
{
	double t, d = 0;
	UINT dim = _dim;
	while (dim--) {
		t = a->x[dim] - b->x[dim];
		d += t * t;
	}
	return d;
}

/**
 * constructor : makes the tree if a valid array of nodes are given
 */
kdtree::kdtree(struct kd_node_t *& t, UINT len, UINT dim) : _dim(dim)
{
    _root = nullptr;
    if ( len > 0)
    {
        _root = make_tree ( t, len, 0); 
    }
}

kdtree::~kdtree()
{
    // no private resources needs deallocation
}

/**
 * swaps contents of a node
 */
void kdtree::swap_node(struct kd_node_t *x, struct kd_node_t * y)
{
    double tmp[MAX_DIM];
    memcpy(tmp,  x->x, sizeof(tmp));
	memcpy(x->x, y->x, sizeof(tmp));
	memcpy(y->x, tmp,  sizeof(tmp));
}



// REFERENCE: http://rosettacode.org/wiki/K-d_tree
