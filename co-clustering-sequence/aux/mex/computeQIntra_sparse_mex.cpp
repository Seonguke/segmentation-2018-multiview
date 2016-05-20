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

MatlabMultiArray<double> compute_color_similarities(MatlabMultiArray<double> H, int length_descriptors, int num_regions, double color_variance);


void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[] )
{
    if(nrhs != 3)
        mexErrMsgTxt("There should be 3 inputs (image, partition and idx_neighbors)");
    
    double *img, *partition, *idx_neighbors_min, *idx_neighbors_max;
    int m_I, n_I, z_I, n_reg;
    mxArray *idx_min_pr, *idx_max_pr;
    

    /* Input parammeters
     * prhs[0] -> img               (double(matrix))
     * prhs[1] -> partittion        (double(matrix))
     * prhs[2] -> idx_neighbors     (struct)
     */
    
    img = mxGetPr(prhs[0]);
    partition = mxGetPr(prhs[1]);
    //idx_neighbors = mxGetPr(prhs[2]);
    
    m_I = mxGetM(prhs[0]);
    n_I = mxGetN(prhs[0])/3;    // We divide by 3 as it is an RGB matrix
    z_I = m_I * n_I;

    
    // Number of regions
    n_reg = 0;
    
    for(int i = 0; i < z_I; i++){
        
        if(partition[i] > n_reg)
            n_reg = partition[i];

    }
    
    
    // Set RGB channel matrices
    mxArray *mat_i_R_pr = mxCreateDoubleMatrix(m_I,n_I,mxREAL);
    mxArray *mat_i_G_pr = mxCreateDoubleMatrix(m_I,n_I,mxREAL);
    mxArray *mat_i_B_pr = mxCreateDoubleMatrix(m_I,n_I,mxREAL);
    
    double *mat_i_R = mxGetPr(mat_i_R_pr);
    double *mat_i_G = mxGetPr(mat_i_G_pr);
    double *mat_i_B = mxGetPr(mat_i_B_pr);
    
    for(int i = 0; i < z_I; i++){
        
        mat_i_R[i] = img[i];
        mat_i_G[i] = img[i+z_I];
        mat_i_B[i] = img[i+2*z_I];
        
    }
    

    int n_bins_rgb = 10;
    
    //Linspace
    mxArray *b_rgb_pr = mxCreateDoubleMatrix(n_bins_rgb+1,1,mxREAL);
    double *b_rgb = mxGetPr(b_rgb_pr);
    double beg = 0;
    double end = 255;
    double step = (end-beg)/n_bins_rgb;
    
    for(int i=0; i<n_bins_rgb+1; i++)
        b_rgb[i] = beg + step*i;
    

    
    /* Output parammeters
     * plhs[0] -> Q     (double(matrix))
     */
    
    mxArray *H_pr = mxCreateDoubleMatrix(3*n_bins_rgb,n_reg,mxREAL);
    MatlabMultiArray<double> H(H_pr);
    
//     plhs[0] = mxCreateDoubleMatrix(3*n_bins_rgb,n_reg,mxREAL);
//     MatlabMultiArray<double> H(plhs[0]);
    
    int n_bin;

    for(int i = 0; i < z_I; i++){
        
        for(int j = 0; j < n_bins_rgb; j++)
        {
            
            if((b_rgb[j] <= mat_i_R[i]) && (mat_i_R[i] < b_rgb[j+1]))
                
                H[j][partition[i]-1]++ ;
            
            if((b_rgb[j] <= mat_i_G[i]) && (mat_i_G[i] < b_rgb[j+1]))
                
                H[j+n_bins_rgb][partition[i]-1]++ ;
            
            if((b_rgb[j] <= mat_i_B[i]) && (mat_i_B[i] < b_rgb[j+1]))
                
                H[j+2*n_bins_rgb][partition[i]-1]++ ;

        }
         
        if(mat_i_R[i] == b_rgb[n_bins_rgb])
            
            H[n_bins_rgb-1][partition[i]-1]++ ;
        
        if(mat_i_G[i] == b_rgb[n_bins_rgb])
            
            H[n_bins_rgb + n_bins_rgb-1][partition[i]-1]++ ;
        
        if(mat_i_B[i] == b_rgb[n_bins_rgb])
            
            H[2*n_bins_rgb + n_bins_rgb-1][partition[i]-1]++ ;
            
        
    }
    
    mxDestroyArray(mat_i_R_pr);
    mxDestroyArray(mat_i_G_pr);
    mxDestroyArray(mat_i_B_pr);
    
    int norm = 0;
    
    for(int jj = 0; jj < n_reg; jj++){
     
        for(int i = 0; i < n_bins_rgb; i++)
            norm = norm + H[i][jj];
        
        
        for(int i = 0; i < 3*n_bins_rgb; i++)
            H[i][jj] = H[i][jj]/norm;
        
        norm = 0;
        
    }
    
    
    double color_variance = 0.4;
    mxArray *u_pr = mxCreateDoubleMatrix(n_reg,n_reg,mxREAL);
    MatlabMultiArray<double> u(u_pr);
    
    idx_min_pr = mxGetField(prhs[2],0,"matrix_min");
    idx_max_pr = mxGetField(prhs[2],0,"matrix_max");
    
    MatlabMultiArray<double> idx_min(idx_min_pr);
    MatlabMultiArray<double> idx_max(idx_max_pr);
    
    int m_idx = mxGetM(idx_max_pr);
    int n_idx = mxGetN(idx_max_pr);
     
    mxArray *v_pr = mxCreateDoubleMatrix(n_reg,n_reg,mxREAL);
    MatlabMultiArray<double> v(v_pr);

    u = compute_color_similarities(H,3*n_bins_rgb,n_reg,color_variance);
    
    int k,l;    
    
    for(int mm = 0; mm < m_idx; mm++){
        for(int nn = 0; nn < n_idx; nn++){
            
            if(idx_max[mm][nn] !=0){
                
                k = idx_max[mm][nn]-1;
                l = idx_min[mm][nn]-1;
                
                v[k][l] = v[k][l]+1;
                v[l][k] = v[k][l];
                
                v[l][l] = v[l][l]+1;
                v[k][k] = v[k][k]+1;
                
            }
        }
    }

    
    double lambda = 0.2;
    int n_diff_zero;
    
    mxArray *Q_pr = mxCreateDoubleMatrix(n_reg,n_reg,mxREAL);
    MatlabMultiArray<double> Q(Q_pr);
    
    for(int i = 0; i < n_reg; i++){
        for(int j = 0; j < n_reg; j++){
            Q[i][j] = lambda * (v[i][j] * u[i][j]);
//             Q[i][j] = lambda * u[i][j];
            
            if(Q[i][j] != 0)
                n_diff_zero++;
        }
    }
    
    plhs[0] = mxCreateSparse(n_reg,n_reg,n_diff_zero,mxREAL); 
    
    mwIndex *irs, *jcs;
    double *sr;
    
    sr = mxGetPr(plhs[0]);
    irs = mxGetIr(plhs[0]);
    jcs = mxGetJc(plhs[0]); 
    
    int idx, accB;
    idx = 0;
    accB = 0;
    
    jcs[idx] = 0;
    idx ++;
    
    for(int nn = 0; nn < n_reg; nn++){
        for(int mm = 0; mm < n_reg; mm++){
            
            if(Q[mm][nn] != 0){
                
                sr[accB] = Q[mm][nn];
                irs[accB] = mm;
                accB++;
            }  
        }
        jcs[idx] = accB;
        idx ++;
    }
    
    

    mxDestroyArray(b_rgb_pr);
    mxDestroyArray(H_pr);
    mxDestroyArray(u_pr);
    mxDestroyArray(v_pr);
    mxDestroyArray(Q_pr);
    
    
    
}



MatlabMultiArray<double> compute_color_similarities(MatlabMultiArray<double> H, int length_descriptors, int num_regions, double color_variance){
       
   
    /* Output allocation */    
    mxArray *u_pr = mxCreateDoubleMatrix(num_regions,num_regions,mxREAL);
    MatlabMultiArray<double> u(u_pr);
    
    double exponent = 0;
    
    for(int k=0; k<num_regions; k++)
    {
        for(int l=0; l<num_regions; l++)
        {
            exponent = 0;
            for(int i=0; i<length_descriptors; i++)
            {
                exponent = exponent + (H[i][k]-H[i][l])*(H[i][k]-H[i][l]);
            }
            exponent = -exponent/(color_variance*color_variance);
            u[k][l] = std::exp(exponent);    
        }
    }
    
    return u;
    
}