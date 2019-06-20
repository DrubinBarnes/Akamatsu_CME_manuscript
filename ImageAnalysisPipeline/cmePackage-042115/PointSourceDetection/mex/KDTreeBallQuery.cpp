 /* [idx, dist] = KDTreeBallQuery(inPts,queryPts,radii);
 *
 * (c) Sylvain Berlemont, 2011 (last modified Jul 14, 2011)
 *
 * Compilation:
 * Mac/Linux: mex -I.  -I../../mex/include/c++ KDTreeBallQuery.cpp
 * Windows: mex COMPFLAGS="$COMPFLAGS /TP /MT" -I"." -I"..\..\mex\include\c++" -output KDTreeBallQuery KDTreeBallQuery.cpp
 */


#include <mex.h>

#include <list>
#include <map>
#include <vector>

#include <vector.hpp>
#include <KDTree.hpp>

template <unsigned K>
static void dispatch(int n, int m, double *x_ptr, double *c_ptr, double *d_ptr, int nlhs, mxArray *plhs[])
{
  // Read parameters
  typename KDTree<K, double>::points_type X;
  
  typename KDTree<K, double>::point_type v;
  
  for (int i = 0; i < n; ++i)
    {
      for (unsigned k = 0; k < K; ++k)
	v[k] = x_ptr[i + (n * k)];
      X.push_back(v);
    }

  typename KDTree<K, double>::points_type C;

  std::vector<double> R;

  for (int i = 0; i < m; ++i)
    {
      for (unsigned k = 0; k < K; ++k)
	v[k] = c_ptr[i + (m * k)];
      C.push_back(v);
      R.push_back(d_ptr[i]);
    }

  // Build kd-tree
  KDTree<K, double> kdtree(X);

  // Compute queries
  std::list<typename KDTree<K, double>::set_type > res_list;
  typename KDTree<K, double>::set_type res;
	
  for (int i = 0; i < m; ++i)
    {
      kdtree.ball_query(C[i], R[i], res);

      res_list.push_back(res);
			
      res.clear();
    }

  // Write output
  if (nlhs > 0)
    {
      const mwSize size[2] = {res_list.size(), 1};
      plhs[0] = mxCreateCellArray(2, size);
			
      int cnt = 0;
      for (typename std::list<typename KDTree<K, double>::set_type>::const_iterator it = res_list.begin(); it != res_list.end(); ++it, ++cnt)
	{
	  if (it->size())
	    {
	      mxArray * pArray = mxCreateDoubleMatrix(it->size(), 1, mxREAL);
	      double * p = mxGetPr(pArray);
				
	      for (typename KDTree<K, double>::set_type::const_iterator it2 = it->begin(); it2 != it->end(); ++it2)
		*p++ = it2->second + 1;
				
	      mxSetCell(plhs[0], cnt, pArray);
	    }
	}
    }

  if (nlhs > 1)
    {
      const mwSize size[2] = {res_list.size(), 1};
      plhs[1] = mxCreateCellArray(2, size);
			
      int cnt = 0;
      for (typename std::list<typename KDTree<K, double>::set_type>::const_iterator it = res_list.begin(); it != res_list.end(); ++it, ++cnt)
	{
	  if (it->size())
	    {
	      mxArray * pArray = mxCreateDoubleMatrix(it->size(), 1, mxREAL);
	      double * p = mxGetPr(pArray);
				
	      for (typename KDTree<K, double>::set_type::const_iterator it2 = it->begin(); it2 != it->end(); ++it2)
		*p++ = it2->first;
				
	      mxSetCell(plhs[1], cnt, pArray);
	    }
	}
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  // Check input/output parameter

  if (nrhs != 3)
    mexErrMsgTxt("Three input arguments required.");

  if (nlhs > 2)
    mexErrMsgTxt("Too many output arguments.");

  size_t k = mxGetN(prhs[0]);
  size_t n = mxGetM(prhs[0]);
  size_t m = mxGetM(prhs[1]);

  if (k != mxGetN(prhs[1]))
    mexErrMsgTxt("X and C must have the same number of columns.");

  if (m != mxGetM(prhs[2]) && mxGetM(prhs[2])!=1)
    mexErrMsgTxt("C and R must have the same number of rows or R must be a scalar.");

  double *x_ptr = mxGetPr(prhs[0]);
  double *c_ptr = mxGetPr(prhs[1]);
  double *d_ptr;
  if (m == mxGetM(prhs[2]))
    d_ptr = mxGetPr(prhs[2]);
  else
    {d_ptr = new double[m];
    for(int i = 0; i < m; ++i) {d_ptr[i] = mxGetScalar(prhs[2]);}
    }

  switch (k)
    {
    case 1: dispatch<1>(n, m, x_ptr, c_ptr, d_ptr, nlhs, plhs); break;
    case 2: dispatch<2>(n, m, x_ptr, c_ptr, d_ptr, nlhs, plhs); break;
    case 3: dispatch<3>(n, m, x_ptr, c_ptr, d_ptr, nlhs, plhs); break;
    default: mexErrMsgTxt("Dimension not implemented.");
    }  
  
  if (m != mxGetM(prhs[2]))
    delete[] d_ptr;
  
}
