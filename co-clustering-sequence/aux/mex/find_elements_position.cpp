#include "mex.h"
#include "matlab_multiarray.hpp"
#include <iostream>
#include <list>
#include <set>
#include <algorithm>

void mexFunction( int nlhs, mxArray *plhs[], 
        		  int nrhs, const mxArray*prhs[] )
{
    double elements;
    if(nrhs<1)
        mexErrMsgTxt("There should be at least 1 input parameter");
    
    /* Input parameters */
    ConstMatlabMultiArray<double> labeled_elements(prhs[0]);  
    elements = mxGetScalar(prhs[1]);  
    
    std::size_t size1   = labeled_elements.shape()[0];
    std::size_t size2   = labeled_elements.shape()[1];
    
    /* Output allocation */    
    plhs[0] = mxCreateDoubleMatrix(elements,2,mxREAL);
    MatlabMultiArray<double> elements_positions(plhs[0]);
    
    for(int i=0; i<size1; i++)
    {
        for(int j=0; j<size2; j++)
        {
            if(labeled_elements[i][j]!=0)
            {
                elements_positions[labeled_elements[i][j]-1][0] = i+1;
                elements_positions[labeled_elements[i][j]-1][1] = j+1;
            }
        }
    }
}