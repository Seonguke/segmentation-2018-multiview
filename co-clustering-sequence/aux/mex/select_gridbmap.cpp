#include "mex.h"
#include "matlab_multiarray.hpp"
#include <iostream>
#include <list>
#include <set>
#include <algorithm>

void mexFunction( int nlhs, mxArray *plhs[], 
        		  int nrhs, const mxArray*prhs[] )
{
    if(nrhs<1)
        mexErrMsgTxt("There should be at least 1 input parameter");
    
    /* Input parameters */
    ConstMatlabMultiArray<double> gridbmap(prhs[0]);
    ConstMatlabMultiArray<double> idx_neighbors_matrix_min(prhs[1]);
    ConstMatlabMultiArray<double> idx_neighbors_matrix_max(prhs[2]);
    ConstMatlabMultiArray<double> v(prhs[3]);
    
    std::size_t size_gridbmap1   = gridbmap.shape()[0];
    std::size_t size_gridbmap2   = gridbmap.shape()[1];
    
    std::size_t idx_neighbors1   = idx_neighbors_matrix_min.shape()[0];
    std::size_t idx_neighbors2   = idx_neighbors_matrix_min.shape()[1];
    
    /* Output allocation */    
    plhs[0] = mxCreateDoubleMatrix(size_gridbmap1,size_gridbmap2,mxREAL);
    MatlabMultiArray<double> gridbmap_out(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(idx_neighbors1,idx_neighbors2,mxREAL);
    MatlabMultiArray<double> idx_neighbors_out_matrix_min(plhs[1]);
    
    plhs[2] = mxCreateDoubleMatrix(idx_neighbors1,idx_neighbors2,mxREAL);
    MatlabMultiArray<double> idx_neighbors_out_matrix_max(plhs[2]);

    plhs[3] = mxCreateDoubleMatrix(size_gridbmap1,size_gridbmap2,mxREAL);
    MatlabMultiArray<double> angles(plhs[3]);
    
    int element = 1, i2, j2;
    for(int i=0; i<size_gridbmap1; i++)
    {
        for(int j=0; j<size_gridbmap2; j++)
        {
            if(gridbmap[i][j]!=0)
            {
                if((i+1) % 2 == 1)
                {
                    i2 = (i+1)/2;
                    j2 = j/2;
                }
                else{
                    i2 = i/2;
                    j2 = (j+1)/2;
                }
                if(angles[i][j]!=0)
                {
                    angles[i][j] = v[i2][j2];
                }
                else
                {
                    if((i+1) % 2 == 1)
                    {
                        angles[i][j]=5;
                    }
                    else{
                        angles[i][j]=1;
                    }
                }
                
//                 if(angles[i][j] != 0)
//                 {
                    gridbmap_out[i][j] = 1;
                    idx_neighbors_out_matrix_min[i][j] = idx_neighbors_matrix_min[i][j];
                    idx_neighbors_out_matrix_max[i][j] = idx_neighbors_matrix_max[i][j];
//                 }
                
            }
        }
    }
}