#include "mex.h"
#include "matlab_multiarray.hpp"
#include <iostream>
#include <list>
#include <set>
#include <algorithm>

void mexFunction( int nlhs, mxArray *plhs[], 
        		  int nrhs, const mxArray*prhs[] )
{
    double variables;
    if(nrhs<1)
        mexErrMsgTxt("There should be at least 1 input parameter");
    
    /* Input parameters */
    ConstMatlabMultiArray<double> RAG(prhs[0]);
    variables = mxGetScalar(prhs[1]);
    
    std::size_t RAG_x = RAG.shape()[0];
    std::size_t RAG_y = RAG.shape()[1];

    /* Output allocation */    
    plhs[0] = mxCreateDoubleMatrix(2, variables, mxREAL);
    MatlabMultiArray<double> regions_adjacency(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(RAG_x, RAG_y, mxREAL);
    MatlabMultiArray<double> regions_variables(plhs[1]);
    regions_variables[0][1] = 1;
    std::cout << "RAG_x:" << RAG_x << std::endl;
    
    double count=0;
    for(int ii=0; ii<RAG_x; ii++)
    {
        for(int jj=0; jj<RAG_y; jj++)
        {
//             std::cout << "RAG[ii][jj]:" << RAG[ii][jj] << std::endl;
            if(RAG[ii][jj] == 1)
            {
                regions_adjacency[0][count] = ii+1;
                regions_adjacency[1][count] = jj+1;
                
                regions_variables[ii][jj] = count+1;
//                 regions_variables[jj][ii] = count+1;
                
                count++;
            }
        }
    }
    
}