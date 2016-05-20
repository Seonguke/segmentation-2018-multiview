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
    
    
    /* Input parameters
     * prhs[0] -> idx_neighbors      (double(matrix))
     * prhs[1] -> ucm                (double(matrix))
     */
    
    double *idx_neighbors_min, *idx_neighbors_max, *ucm;
    mxArray *idx_neighbors_min_pr, *idx_neighbors_max_pr;
    int N_reg;
    
    idx_neighbors_min_pr = mxGetField(prhs[0],0,"matrix_min");
    idx_neighbors_max_pr = mxGetField(prhs[0],0,"matrix_max");
    
    idx_neighbors_min = mxGetPr(idx_neighbors_min_pr);
    idx_neighbors_max = mxGetPr(idx_neighbors_max_pr);
    
    int m_grid = mxGetM(idx_neighbors_min_pr);
    int n_grid = mxGetN(idx_neighbors_min_pr);

    ucm = mxGetPr(prhs[1]);

    N_reg = 0;
    
    for(int j = 0; j < n_grid; j++){
        for(int i =0; i< m_grid; i++){
            if(idx_neighbors_max[i + j*m_grid]>N_reg){
                N_reg = idx_neighbors_max[i + j*m_grid];
            }
        }
    }

    
    /* Output allocation
     * plhs[0] -> RAG_Intra         (double(matrix))
     * plhs[1] -> RAG_Inter         (double(matrix))
     */
    
    
    plhs[0] = mxCreateDoubleMatrix(N_reg,N_reg,mxREAL);
    MatlabMultiArray<double> weights_QIntra(plhs[0]);
    
    

    // INTER
    for(int j = 0; j < n_grid; j++){
        for(int i =0; i< m_grid; i++){
            if(idx_neighbors_max[i + j*m_grid]!=0){
                int idx_reg_1 = idx_neighbors_min[i + j*m_grid];
                int idx_reg_2 = idx_neighbors_max[i + j*m_grid];
                //weights_QIntra[idx_reg_1-1][idx_reg_2-1] = 1-ucm[i + j*m_grid];
                //weights_QIntra[idx_reg_2-1][idx_reg_1-1] = 1-ucm[i + j*m_grid];
                weights_QIntra[idx_reg_1-1][idx_reg_2-1] = exp(-1*ucm[i + j*m_grid]/0.2);
                weights_QIntra[idx_reg_2-1][idx_reg_1-1] = exp(-1*ucm[i + j*m_grid]/0.2);
            }
        }
    }
    
    
}