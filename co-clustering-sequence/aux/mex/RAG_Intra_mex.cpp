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

using namespace std;    


void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[] )
{
    if(nrhs != 2)
        mexErrMsgTxt("There should be 2 inputs (idx_neighbor and partition)");
    
    
    
    double  *partition1, *idx_neighbors_min, *idx_neighbors_max;
    mxArray *idx_neighbors_min_pr, *idx_neighbors_max_pr;
    int N_reg1;
    
    
    /* Input parameters
     * prhs[0] -> idx_neighbor   (structure)
     * prhs[1] -> partition1     (double(matrix))
     */
    
    
    partition1 = mxGetPr(prhs[1]);
    
    int m = mxGetM(prhs[1]);
    int n = mxGetN(prhs[1]);
    
    
//     mxArray *partition1_data_pr = mxCreateDoubleMatrix(m,n,mxREAL);
//     MatlabMultiArray<double> partition1_data(partition1_data_pr);
   
    
    N_reg1 = 0;
    
    for(int i = 0; i < m*n; i++){
   
        if(partition1[i] > N_reg1)
            N_reg1 = partition1[i];

    }
    
    
    /* Output allocation
     * plhs[0] -> RAG_Intra         (double(matrix))
     */
    
    plhs[0] = mxCreateDoubleMatrix(N_reg1,N_reg1,mxREAL);
    MatlabMultiArray<double> RAG_Intra(plhs[0]);
    
    
    // INTRA
    
    idx_neighbors_min_pr = mxGetField(prhs[0],0,"matrix_min");
    idx_neighbors_max_pr = mxGetField(prhs[0],0,"matrix_max");
    
    idx_neighbors_min = mxGetPr(idx_neighbors_min_pr);
    idx_neighbors_max = mxGetPr(idx_neighbors_max_pr);
    
    int m_grid = mxGetM(idx_neighbors_min_pr);
    int n_grid = mxGetN(idx_neighbors_min_pr);
    
    
    for(int j = 0; j < n_grid; j++){
        for(int i =0; i< m_grid; i++){
            
            if((idx_neighbors_min[i + j*m_grid] != 0) && (idx_neighbors_max[i + j*m_grid] != 0)){
                
                RAG_Intra[idx_neighbors_min[i + j*m_grid]-1][idx_neighbors_max[i + j*m_grid]-1] = 1;
                RAG_Intra[idx_neighbors_max[i + j*m_grid]-1][idx_neighbors_min[i + j*m_grid]-1] = 1;
                
            } 
        }
    }
    

    
}
