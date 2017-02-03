/*
------------------------------------------------------------------------
  Copyright (C)
  Universitat Politecnica de Catalunya BarcelonaTech (UPC) - Spain

  Andreu Girbau <andreu.girbau@alu-etsetb.upc.edu>
  September 2014
------------------------------------------------------------------------
 */


#include "mex.h"
#include "matlab_multiarray.hpp"
#include <iostream>
#include <list>
#include <set>
#include <algorithm>
#include "matrix.h"
#include "string.h"

using namespace std;    // For std::vector



void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[] )
{
    if(nrhs != 2)
        mexErrMsgTxt("There should be 2 inputs (partition1 and partition2)");
       
    
    double *partition1, *partition2;
    int N_reg1, N_reg2;
    
    
    /* Input parameters
     * prhs[0] -> partition      (double(matrix))
     * prhs[1] -> partition2     (double(matrix))
     */
    
    
    partition1 = mxGetPr(prhs[0]);
    partition2 = mxGetPr(prhs[1]);
    
    int m = mxGetM(prhs[0]);
    int n = mxGetN(prhs[0]);
    
    
    mxArray *partition1_data_pr = mxCreateDoubleMatrix(m,n,mxREAL);
    MatlabMultiArray<double> partition1_data(partition1_data_pr);
    
    mxArray *partition2_data_pr = mxCreateDoubleMatrix(m,n,mxREAL);
    MatlabMultiArray<double> partition2_data(partition2_data_pr);

    
    N_reg1 = 0;
    N_reg2 = 0;
    
    
    for(int i = 0; i < m*n; i++){
        
        if(partition1[i] > N_reg1)
            N_reg1 = partition1[i];

        if(partition2[i] > N_reg2)
            N_reg2 = partition2[i];  
        
    }

    
    /* Output allocation
     * plhs[0] -> RAG_Intra         (double(matrix))
     * plhs[1] -> RAG_Inter         (double(matrix))
     */
    
    
    plhs[0] = mxCreateDoubleMatrix(N_reg1,N_reg2,mxREAL);
    MatlabMultiArray<double> RAG_Inter(plhs[0]);
    

    // INTER
    
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
             RAG_Inter[partition1[i + j*m]-1][partition2[i + j*m]-1] = 1;
             
        }
    }
    
    
}
