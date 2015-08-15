#include "kdtree.h"

void kdtree::nearest(struct kd_node_t *root, struct kd_node_t *nd, int i,
		struct kd_node_t **best, double *best_dist)
{
	double d, dx, dx2;
 
	if (!root) return;
	d = dist(root, nd);
	dx = root->x[i] - nd->x[i];
	dx2 = dx * dx;

	if (!*best || d < *best_dist) 
	{
		*best_dist = d;
		*best = root;
	}
 
	/* if chance of exact match is high */
	if (!*best_dist) return;
 
	if (++i >= _dim) i = 0;
 
	nearest(dx > 0 ? root->left : root->right, nd, i, best, best_dist);
	if (dx2 >= *best_dist) return;
	nearest(dx > 0 ? root->right : root->left, nd, i, best, best_dist);
}

struct kd_node_t* 
kdtree::make_tree(struct kd_node_t *t, int len, int i)
{
	struct kd_node_t *n;
 
	if (!len) return 0;
 
	if ((n = find_median(t, t + len, i))) {
		i = (i + 1) % _dim;
		n->left  = make_tree(t, n - t, i)
		n->right = make_tree(n + 1, t + len - (n + 1), i);
	}
	return n;
}

struct kd_node_t*
find_median(struct kd_node_t *start, struct kd_node_t *end, int idx)
{
	if (end <= start) return NULL;
	if (end == start + 1)
		return start;
 
	struct kd_node_t *p, *store, *md = start + (end - start) / 2;
	double pivot;
	while (1) {
		pivot = md->x[idx];
 
		swap(md, end - 1);
		for (store = p = start; p < end; p++) {
			if (p->x[idx] < pivot) {
				if (p != store)
					swap(p, store);
				store++;
			}
		}
		swap(store, end - 1);
 
		/* median has duplicate values */
		if (store->x[idx] == md->x[idx])
			return md;
 
		if (store > md) end = store;
		else        start = store;
	}
}

double
kdtree::dist(struct kd_node_t *a, struct kd_node_t *b)
{
	double t, d = 0;
	UINT dim = _dim;
	while (dim--) {
		t = a->x[dim] - b->x[dim];
		d += t * t;
	}
	return d;
}

kdtree::kdtree(UINT dim) : _dim(dim)
{
	_root = make_tree ()	
}

// REFERENCE: http://rosettacode.org/wiki/K-d_tree