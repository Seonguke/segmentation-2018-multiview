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
#include "math.h"
#include <complex>

using namespace std;    

complex<double> mycomplex (0.0,1.0);
double PI_trans = 3.1416/180;


void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[] )
{
    if(nrhs != 6)
        mexErrMsgTxt("There should be 6 inputs (labeled_descriptors1, labeled_descriptors2, labeled_elements1, labeled_elements2, Sigma and window)");
    
    
    
    //double  *labeled_descriptors1, *labeled_descriptors2, *labeled_elements1, *labeled_elements2, *sigma;
    int window;
    int m_le, n_le, z_le, m_ld, n_ld1, n_ld2, accW_value;
    double accum;
    vector<double> value, accW;
    vector<int> m_W;
//     double *accW;
    
    
    /* Input parammeters
     * prhs[0] -> labeled_descriptors1    (double(matrix))
     * prhs[1] -> labeled_descriptors2    (double(matrix))
     * prhs[2] -> labeled_elements1       (double(matrix))
     * prhs[3] -> labeled_elements2       (double(matrix))
     * prhs[4] -> sigma                   (double(matrix))
     * prhs[5] -> window                  (int)
     */
    
//     labeled_descriptors1 = mxGetPr(prhs[0]);
//     labeled_descriptors2 = mxGetPr(prhs[1]);
//     labeled_elements1 = mxGetPr(prhs[2]);
//     labeled_elements2 = mxGetPr(prhs[3]); 
//     sigma = mxGetPr(prhs[4]);
     window = mxGetScalar(prhs[5]);
    
    ConstMatlabMultiArray<double> labeled_descriptors1(prhs[0]);
    ConstMatlabMultiArray<double> labeled_descriptors2(prhs[1]);
    ConstMatlabMultiArray<double> labeled_elements1(prhs[2]);
    ConstMatlabMultiArray<double> labeled_elements2(prhs[3]);
    ConstMatlabMultiArray<double> Sigma(prhs[4]);

    m_le = mxGetM(prhs[2]);
    n_le = mxGetN(prhs[2]);
    z_le = m_le * n_le;

    m_ld = mxGetM(prhs[1]);
    
    n_ld1 = mxGetN(prhs[0]);
    n_ld2 = mxGetN(prhs[1]);
    

    int mi,ni,mj,nj,d, n_diff_zero;
    double dist;
    
    accW_value = 0;
    n_diff_zero = 0;
    
    //Buscarem des de labeled_2 a labeled_1 per l'sparse matrix en HORITZONTAL (per això està definit en matriu, per facilitat)
    
    for(int ii = 0; ii < m_le; ii++){
        for(int jj = 0; jj < n_le; jj++){
            
            if(labeled_elements2[ii][jj] != 0){
                // Buscar elements circundants a labeled_elements2[ii]
                
                for(int mm = -window; mm <= window; mm++){
                    for(int nn = -window; nn <= window; nn++){
                        
                        // Comprovació per veure si l'element que busquem està a l'interior
                        if((jj+nn >= 0) && (jj+nn < n_le) && (ii+mm >= 0) && (ii+mm < m_le)){ 
                            
                            if(labeled_elements1[ii+mm][jj+nn] != 0){
                                
                                accum = 0;
                                
                                mi = ii;
                                ni = jj;
                                mj = ii+mm;
                                nj = jj+nn;
                                
                                dist = sqrt(((double(mi)-double(mj))*(double(mi)-double(mj))+(double(ni)-double(nj))*(double(ni)-double(nj))));
                                
                                if(dist <= window){
                                    
                                    for(int kk = 0; kk < m_ld; kk++){
                                        
                                        accum = accum + (Sigma[kk][kk] *
                                                (labeled_descriptors2[kk][(int)labeled_elements2[ii][jj]-1] - labeled_descriptors1[kk][(int)labeled_elements1[ii+mm][jj+nn]-1])*
                                                (labeled_descriptors2[kk][(int)labeled_elements2[ii][jj]-1] - labeled_descriptors1[kk][(int)labeled_elements1[ii+mm][jj+nn]-1]));
                                        
                                    }
                                    
                                    accum = -1*accum;
                                    
                                    value.push_back(exp(accum));
                                    m_W.push_back(labeled_elements1[ii+mm][jj+nn]-1);
 
                                    n_diff_zero++;
                                    
                                }
                            }
                        }
                    }
                }
                
                 accW.push_back(value.size());
                
            }
        }
    }

  
    /* Output parammeters
     * plhs[0] -> labeled_element_descriptors     (sparse(matrix))
     */
    
    plhs[0] = mxCreateSparse(n_ld1,n_ld2,n_diff_zero,mxREAL);
    
    mwIndex *irs, *jcs;
    double *sr;
    
    sr  = mxGetPr(plhs[0]);
    irs = mxGetIr(plhs[0]);
    jcs = mxGetJc(plhs[0]); 
    
    jcs[0] = 0;
      
    
    for(int i = 0; i < n_diff_zero; i++){    

        sr[i] = value[i];
        irs[i] = m_W[i];       

    } 

    
    for(int i = 0; i < n_ld2; i++){
     
        jcs[i+1] = accW[i];
        
    }
    
 
 
}