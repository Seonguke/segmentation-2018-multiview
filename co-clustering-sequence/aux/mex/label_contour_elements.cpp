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
    
    std::size_t size_gridbmap1   = gridbmap.shape()[0];
    std::size_t size_gridbmap2   = gridbmap.shape()[1];
    
    /* Output allocation */    
    plhs[0] = mxCreateDoubleMatrix(size_gridbmap1,size_gridbmap2,mxREAL);
    MatlabMultiArray<double> labeled_elements(plhs[0]);
    
    int element = 1;
    for(int i=0; i<size_gridbmap1; i++)
    {
        for(int j=0; j<size_gridbmap2; j++)
        {
            if(gridbmap[i][j]!=0)
            {
                labeled_elements[i][j] = element;
                element++;
            }
        }
    }
}
