#ifndef KDTREE_HPP
# define KDTREE_HPP

# include <map>
# include <vector>
# include <algorithm>

# include <vector.hpp>

#if defined(_WIN32) || defined(_WIN64)
double log2(double n) { return log(n) / log(2.0); }
#endif

template <unsigned N, typename T>
class KDTree
{
public:
  typedef vector<N, T> point_type;
  typedef typename std::vector<point_type> points_type;
  typedef typename std::vector< unsigned >::iterator idx_iterator_type;
  typedef std::pair<double, unsigned> pair_type;
  typedef std::multimap<double, unsigned> set_type;
	
public:
  KDTree(const points_type & points) : points_(points), tree_(NULL)
  {
    int n = points_.size();

    if (n == 0)
      return;
		
    // n must be less than 2^32 - 1 since indices are stored on signed int
    assert(n <= std::numeric_limits<int>::max());
		
    std::vector<unsigned> indices(n);
    std::generate(indices.begin(), indices.end(), Incr_());
		
    int h = (int) ceil(log2(n)) + 1;
    int nnodes = (1 << (h+1)) - 1;
    tree_ = new int[nnodes];
		
    build_(indices.begin(), indices.end());
  }

  KDTree(const points_type & points, const std::vector<int> & tree)
  {
    int n = points_.size();

    if (n == 0)
      return;

    // n must be less than 2^32 - 1 since indices are stored on signed int
    assert(n <= std::numeric_limits<int>::max());

    int h = (int) ceil(log2(n)) + 1;
    int nnodes = (1 << (h+1)) - 1;

    // input tree must at least has the correct size
    assert(tree.size() == nnodes);

    tree_ = new int[nnodes];
    std::copy(tree.begin(), tree.end(), tree_);
  }

  ~KDTree()
  {
    delete[] tree_;
  }
	
  pair_type closest_point(const point_type & query, unsigned idx = 0, unsigned k = 0) const
  {
    unsigned idx_lt = (idx << 1) + 1;
    unsigned idx_gt = (idx << 1) + 2;
    unsigned k2 = (k + 1) & (N - 1);
		
    if (tree_[idx_lt] != -1)
      {
	double value = points_[tree_[idx]][k];
						
	// top-down
	pair_type res = query[k] < value ?
	  closest_point(query, idx_lt, k2) :
	  closest_point(query, idx_gt, k2);
		
	// bottom-up
	if (query[k] < value)
	  {
	    // The closest point has been found in the lt_ child. See
	    // whether there isn't any closer point in the gt_ child
	    // within min_dist from the closest point.
	    if (value - query[k] < res.first)
	      {
		pair_type res2 = closest_point(query, idx_gt, k2);
					
		if (res2.first < res.first)
		  res = res2;
	      }
	  }
	else
	  {
	    // The closest point has been found in the gt_ child. See
	    // whether there isn't any closer point in the lt_ child
	    // within min_dist from the closest point.
	    if (query[k] - value < res.first)
	      {
		pair_type res2 = closest_point(query, idx_lt, k2);
					
		if (res2.first < res.first)
		  res = res2;
	      }
	  }
			
	return res;
      }
    else
      {
	// idx is a leaf
	return std::make_pair(dist(query, points_[tree_[idx]]), tree_[idx]);
      }
  }
	
  void ball_query(const point_type & center, double radius, set_type & res, unsigned idx = 0, unsigned k = 0) const
  {
    unsigned idx_lt = (idx << 1) + 1;
    unsigned idx_gt = (idx << 1) + 2;
    unsigned k2 = (k + 1) & (N - 1);

    if (tree_[idx_lt] != -1)
      {
	double value = points_[tree_[idx]][k];
			
	double d = center[k] - value;
			
	if (d > 0)
	  {
	    if (d <= radius)
	      ball_query(center, radius, res, idx_lt, k2);
	    ball_query(center, radius, res, idx_gt, k2);
	  }
	else
	  {
	    if (-d <= radius)
	      ball_query(center, radius, res, idx_gt, k2);
	    ball_query(center, radius, res, idx_lt, k2);
	  }			
      }
    else
      {
	//idx is a leaf
	double d = dist(center, points_[tree_[idx]]);
			
	if (d <= radius)
	  res.insert(std::make_pair(d, tree_[idx]));
      }
  }

  void range_query(const point_type & center, const point_type & hside, set_type & res, unsigned idx = 0, unsigned k = 0) const
  {
    unsigned idx_lt = (idx << 1) + 1;
    unsigned idx_gt = (idx << 1) + 2;
    unsigned k2 = (k + 1) & (N - 1);

    if (tree_[idx_lt] != -1)
      {
	double value = points_[tree_[idx]][k];

	double d = center[k] - value;

	if (d > 0)
	  {
	    if (d <= hside[k])
	      range_query(center, hside, res, idx_lt, k2);
	    range_query(center, hside, res, idx_gt, k2);
	  }
	else
	  {
	    if (-d <= hside[k])
	      range_query(center, hside, res, idx_gt, k2);
	    range_query(center, hside, res, idx_lt, k2);
	  }
      }
    else
      {
	//idx is a leaf
	bool isInside = true;

	for (unsigned i = 0; i < N; ++i)
	  isInside &= fabs(center[i] - points_[tree_[idx]][i]) <= hside[i];

	if (isInside)
	  res.insert(std::make_pair(dist(center, points_[tree_[idx]]), tree_[idx]));
      }
  }
	
private:
    void build_(const idx_iterator_type & first,
                const idx_iterator_type & last,
                unsigned idx = 0,
                unsigned k = 0)
    {
        unsigned idx_lt = (idx << 1) + 1;
        unsigned idx_gt = (idx << 1) + 2;
        
        if (last - first > 1)
        {
            unsigned k2 = (k + 1) & (N - 1);
            
            // find the median value:
            idx_iterator_type median = first + (last - first) / 2;
            nth_element (first, median, last, Cmp_(points_, k));
            tree_[idx] = *median;

            // build children
            build_(first, median, idx_lt, k2);
            build_(median, last, idx_gt, k2);            
        }
        else
        {
            tree_[idx] = *first;
            tree_[idx_lt] = -1;
            tree_[idx_gt] = -1;            
        }
    }
	
private:
  struct Incr_
  {
    Incr_() { current = 0; }
		
    int operator()() { return current++; }

    int current;		
  };

  struct Cmp_
  {
  public:
    Cmp_(const points_type & points, unsigned k) : points_(points), k_(k) {}
		
    bool operator()(unsigned i, unsigned j) const
    {
      return points_[i][k_] < points_[j][k_];
    }
		
  private:
    const points_type & points_;
		
    unsigned k_;
  };
	
private:
  const points_type & points_;
	
  int* tree_;
};

#endif /* KDTREE_HPP */
