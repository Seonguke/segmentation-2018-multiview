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

using namespace std;    // For std::vector

complex<double> mycomplex (0.0,1.0);
double PI_trans = 3.1416/180;
int window = 5;
int window_B = window*window;


//MatlabMultiArray<double> hist_c(double *v, int size_v, double *bin_edges, int n_bins);
double *hist_c(double *v, int size_v, double *bin_edges, int n_bins);
//int norm(MatlabMultiArray<double> hist_array, int n_bins);
int norm(double *hist_array, int n_bins);


void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[] )
{
//     if(nrhs != 3)
//         mexErrMsgTxt("There should be 3 inputs (image, gPb_angles and labled_elements)");
    
    
    
    double  *img, *gPb_angles, *labeled_elements, *labeled_element_descriptors, *matrix_min, *matrix_max;
    double  *mat_i_R, *mat_i_G, *mat_i_B;
    int m_I, n_I, z_I, m_B, n_B, z_B, els, n_bins_rgb, n_bins_grad, norm_rgb, norm_grad;
    
    
    
    /* Input parammeters
     * prhs[0] -> img               (double(matrix))
     * prhs[1] -> gPb_angles        (double(matrix))
     * prhs[2] -> labeled_elements  (double(matrix))
     */
    
    img              = mxGetPr(prhs[0]);
    gPb_angles       = mxGetPr(prhs[1]);
    labeled_elements = mxGetPr(prhs[2]);
    matrix_min       = mxGetPr(prhs[3]);
    matrix_max       = mxGetPr(prhs[4]);
    
    m_I = mxGetM(prhs[0]);
    n_I = mxGetN(prhs[0])/3;    // We divide by 3 as it is an RGB matrix
    z_I = m_I * n_I;
    
    m_B = mxGetM(prhs[1]);
    n_B = mxGetN(prhs[1]);
    z_B = m_B * n_B;        //m_B = 2*m_I+1     n_B = 2*n_I+1
    
    
    // Set RGB channel matrices
    mxArray *mat_i_R_pr_mat = mxCreateDoubleMatrix(m_B,n_B,mxREAL);
    mxArray *mat_i_G_pr_mat = mxCreateDoubleMatrix(m_B,n_B,mxREAL);
    mxArray *mat_i_B_pr_mat = mxCreateDoubleMatrix(m_B,n_B,mxREAL);
    
    MatlabMultiArray<double> mat_i_R_mat(mat_i_R_pr_mat);
    MatlabMultiArray<double> mat_i_G_mat(mat_i_G_pr_mat);
    MatlabMultiArray<double> mat_i_B_mat(mat_i_B_pr_mat);
    
    int idx_image = 0;
    
    for(int n=1; n < n_B; n=n+2){
        for(int m=1; m < m_B; m=m+2){
            
            mat_i_R_mat[m][n] = img[idx_image];
            mat_i_G_mat[m][n] = img[idx_image+z_I];
            mat_i_B_mat[m][n] = img[idx_image+2*z_I];
            
            idx_image++;
            
        }
    }
    
    mxArray *mat_i_R_pr = mxCreateDoubleMatrix(m_B,n_B,mxREAL);
    mxArray *mat_i_G_pr = mxCreateDoubleMatrix(m_B,n_B,mxREAL);
    mxArray *mat_i_B_pr = mxCreateDoubleMatrix(m_B,n_B,mxREAL);
    
    mat_i_R = mxGetPr(mat_i_R_pr);
    mat_i_G = mxGetPr(mat_i_G_pr);
    mat_i_B = mxGetPr(mat_i_B_pr);
    
    idx_image = 0;
    
    for(int n=0; n < n_B; n++){
        for(int m=0; m < m_B; m++){
            
            mat_i_R[idx_image] = mat_i_R_mat[m][n];
            mat_i_G[idx_image] = mat_i_G_mat[m][n];
            mat_i_B[idx_image] = mat_i_B_mat[m][n];
            
            idx_image++;
            
        }
    }
    
    mxDestroyArray(mat_i_R_pr_mat);
    mxDestroyArray(mat_i_G_pr_mat);
    mxDestroyArray(mat_i_B_pr_mat);
    
    
    //Compute number of elements
    
    els = 0;
    
    for(int i=0; i < z_B; i++){
        if(els < labeled_elements[i])
            els = labeled_elements[i];
        
    }
    

    n_bins_rgb = 8;
    
    //Linspace 
    mxArray *b_rgb_pr = mxCreateDoubleMatrix(n_bins_rgb+1,1,mxREAL);
    double *b_rgb = mxGetPr(b_rgb_pr);
    double beg = 0;
    double end = 255;
    double step = (end-beg)/n_bins_rgb;
    
    for(int i=0; i<n_bins_rgb+1; i++)
        b_rgb[i] = beg + step*i;
    
    
    //int window = 5;
    //int window_B = window*window;
    
    mxArray *r_surr_pr = mxCreateDoubleMatrix(2*window_B,1,mxREAL);
    mxArray *g_surr_pr = mxCreateDoubleMatrix(2*window_B,1,mxREAL);
    mxArray *b_surr_pr = mxCreateDoubleMatrix(2*window_B,1,mxREAL);
    
    double *r_surr = mxGetPr(r_surr_pr);
    double *g_surr = mxGetPr(g_surr_pr);
    double *b_surr = mxGetPr(b_surr_pr);
    
    n_bins_grad = 8;
    
    //Linspace
    mxArray *b_grad_pr = mxCreateDoubleMatrix(n_bins_grad+1,1,mxREAL);
    double *b_grad = mxGetPr(b_grad_pr);
    beg = 22.5;
    end = 180;
    step = (end-beg)/n_bins_rgb;
    
    for(int i=0; i<n_bins_grad+1; i++)
        b_grad[i] = beg + step*i;
    
    mxArray *grad_surr_pr = mxCreateDoubleMatrix(window_B,1,mxREAL);
    double *grad_surr = mxGetPr(grad_surr_pr);
    
    
    mxArray *hist_r_pr = mxCreateDoubleMatrix(1,n_bins_rgb,mxREAL);
    mxArray *hist_g_pr = mxCreateDoubleMatrix(1,n_bins_rgb,mxREAL);
    mxArray *hist_b_pr = mxCreateDoubleMatrix(1,n_bins_rgb,mxREAL);
    
    double *hist_r = mxGetPr(hist_r_pr);
    double *hist_g = mxGetPr(hist_g_pr);
    double *hist_b = mxGetPr(hist_b_pr);
    
//     MatlabMultiArray<double> hist_r(hist_r_pr);
//     MatlabMultiArray<double> hist_g(hist_g_pr);
//     MatlabMultiArray<double> hist_b(hist_b_pr);
    
    mxArray *hist_grad_pr = mxCreateDoubleMatrix(1,n_bins_grad,mxREAL);
    double *hist_grad = mxGetPr(hist_grad_pr);
//     MatlabMultiArray<double> hist_grad(hist_grad_pr);
    
    
    /* Output parammeters
     * plhs[0] -> labeled_element_descriptors     (double(matrix))
     * plhs[1] -> 3*n_bins_rgb                    (double(scalar))
     * plhs[2] -> n_bins_grad                     (double(scalar))
     */
    
    int rows_plhs = 3*n_bins_rgb+n_bins_grad;
    
    plhs[0] = mxCreateDoubleMatrix(rows_plhs,els,mxREAL);
    labeled_element_descriptors = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    double *out_n_bins_rgb = mxGetPr(plhs[1]);
    
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    double *out_n_bins_grad = mxGetPr(plhs[2]);
    
    out_n_bins_rgb[0] = 3*n_bins_rgb;
    out_n_bins_grad[0] = n_bins_grad;
    
//     plhs[0] = mxCreateDoubleMatrix(1,n_bins_grad,mxREAL);
//     MatlabMultiArray<double> hist_grad(plhs[0]);
    
//     plhs[0] = mxCreateDoubleMatrix(1,n_bins_rgb,mxREAL);
//     MatlabMultiArray<double> hist_r(plhs[0]);
    
    // La manera de treure l'output serà fent for(i fins a 3*n_bins_rgb+n_bins_grad)
    // labeled_element_descriptors(i+labeled_element[i]*(3*n_bins_rgb+n_bins_grad))
    
    
    int idx = 0;
    int mm_grad = 0;
    int nn_grad = 0;
    
    for(int ii=0; ii<z_B; ii++){
//          int ii = 283043;
        if(labeled_elements[ii] != 0){
            
            idx = 0;
            
            if((((ii+1)%m_B)%2) != 0){
                //Horizontal contours
                
                mm_grad = 0;
                for(int mm=0; mm<2*window; mm=mm+2){
                    nn_grad = 0;
                    for(int nn=0; nn<2*window; nn=nn+2){
                        
                        if((ii-1-m_B*(window-1-nn)-(window-1-mm) > 0) && (ii+1-m_B*(window-1-nn)-(window-1-mm) < z_B)){       // Search inside the matrix
                            
                            r_surr[idx]          = mat_i_R[ii-1-m_B*(window-1-nn)-(window-1-mm)];
                            r_surr[idx+window_B] = mat_i_R[ii+1-m_B*(window-1-nn)-(window-1-mm)];
                            g_surr[idx]          = mat_i_G[ii-1-m_B*(window-1-nn)-(window-1-mm)];
                            g_surr[idx+window_B] = mat_i_G[ii+1-m_B*(window-1-nn)-(window-1-mm)];
                            b_surr[idx]          = mat_i_B[ii-1-m_B*(window-1-nn)-(window-1-mm)];
                            b_surr[idx+window_B] = mat_i_B[ii+1-m_B*(window-1-nn)-(window-1-mm)];
                            
                        }
                        else{
                            
                            r_surr[idx]          = 0;
                            r_surr[idx+window_B] = 0;
                            g_surr[idx]          = 0;
                            g_surr[idx+window_B] = 0;
                            b_surr[idx]          = 0;
                            b_surr[idx+window_B] = 0;
                            
                            
                        }
                        
                        if(((ii-m_B*(((window-1)/2)-nn_grad)-(((window-1)/2)-mm_grad)) > 0)&&((ii-m_B*(((window-1)/2)-nn_grad)-(((window-1)/2)-mm_grad)) < z_B)){
//                             std::cout << "Position horizontal: "<<labeled_elements[ii]<< " " <<labeled_elements[ii-m_B*(((window-1)/2)-nn_grad)-(((window-1)/2)-mm_grad)]<<std::endl;
                            if( (matrix_min[ii-m_B*(((window-1)/2)-nn_grad)-(((window-1)/2)-mm_grad)] == matrix_min[ii]) || (matrix_max[ii-m_B*(((window-1)/2)-nn_grad)-(((window-1)/2)-mm_grad)] == matrix_max[ii]))
                            {
                                grad_surr[idx] = gPb_angles[ii-m_B*(((window-1)/2)-nn_grad)-(((window-1)/2)-mm_grad)];
                            }else
                                grad_surr[idx] = -22.5;
                            
                            // grad_surr(idx) = gPb_angles(ii-m_l*(((window-1)/2)-nn_grad+1)-(((window-1)/2)-mm_grad+1));

                        }
                        else
                            grad_surr[idx] = -22.5;
                        
                        idx++;
                        nn_grad++;
                    }
                    
                    mm_grad++;
                }
                
            }
            else{
                //Vertical contours
                
                mm_grad = 0;
                for(int mm=0; mm<2*window; mm=mm+2){
                    nn_grad = 0;
                    for(int nn=0; nn<2*window; nn=nn+2){
                        
                        if((ii-m_B-m_B*(window-1-nn)-(window-1-mm) > 0) && (ii+m_B-m_B*(window-1-nn)-(window-1-mm) < z_B)){       // Search inside the matrix
                            
                            r_surr[idx]          = mat_i_R[ii-m_B-m_B*(window-1-nn)-(window-1-mm)];
                            r_surr[idx+window_B] = mat_i_R[ii+m_B-m_B*(window-1-nn)-(window-1-mm)];
                            g_surr[idx]          = mat_i_G[ii-m_B-m_B*(window-1-nn)-(window-1-mm)];
                            g_surr[idx+window_B] = mat_i_G[ii+m_B-m_B*(window-1-nn)-(window-1-mm)];
                            b_surr[idx]          = mat_i_B[ii-m_B-m_B*(window-1-nn)-(window-1-mm)];
                            b_surr[idx+window_B] = mat_i_B[ii+m_B-m_B*(window-1-nn)-(window-1-mm)];
                            
                            
                        }
                        else{
                            
                            r_surr[idx]          = 0;
                            r_surr[idx+window_B] = 0;
                            g_surr[idx]          = 0;
                            g_surr[idx+window_B] = 0;
                            b_surr[idx]          = 0;
                            b_surr[idx+window_B] = 0;
                            
                            
                        }
                        
                        if(((ii-m_B*(((window-1)/2)-nn_grad)-(((window-1)/2)-mm_grad)) > 0)&&((ii-m_B*(((window-1)/2)-nn_grad)-(((window-1)/2)-mm_grad)) < z_B)){

                            if( (matrix_min[ii-m_B*(((window-1)/2)-nn_grad)-(((window-1)/2)-mm_grad)] == matrix_min[ii]) || (matrix_max[ii-m_B*(((window-1)/2)-nn_grad)-(((window-1)/2)-mm_grad)] == matrix_max[ii]))
                            {
                                grad_surr[idx] = gPb_angles[ii-m_B*(((window-1)/2)-nn_grad)-(((window-1)/2)-mm_grad)];
                            }else
                                grad_surr[idx] = -22.5;
                            

                        }
                        else
                            grad_surr[idx] = -22.5;
                        
                        idx++;
                        nn_grad++;
                    }
                    
                    mm_grad++;
                }
                
                
            }
            
            
            hist_r = hist_c(r_surr,2*window_B,b_rgb,n_bins_rgb+1);
            hist_g = hist_c(g_surr,2*window_B,b_rgb,n_bins_rgb+1);
            hist_b = hist_c(b_surr,2*window_B,b_rgb,n_bins_rgb+1);
            hist_grad = hist_c(grad_surr,window_B,b_grad,n_bins_grad+1);
            
            // Normalization
            norm_rgb = norm(hist_r,n_bins_rgb);
            norm_grad = norm(hist_grad,n_bins_grad);
            
            for(int zz=0; zz<n_bins_rgb; zz++){
                
                                                
                labeled_element_descriptors[zz+((int)labeled_elements[ii]-1)*rows_plhs] = hist_r[zz]/norm_rgb;
                labeled_element_descriptors[zz+n_bins_rgb+((int)labeled_elements[ii]-1)*rows_plhs] = hist_g[zz]/norm_rgb;
                labeled_element_descriptors[zz+2*n_bins_rgb+((int)labeled_elements[ii]-1)*rows_plhs] = hist_b[zz]/norm_rgb;

//                 
//                 labeled_element_descriptors[zz+((int)labeled_elements[ii]-1)*rows_plhs] = hist_r[0][zz]/norm_rgb;
//                 labeled_element_descriptors[zz+n_bins_rgb+((int)labeled_elements[ii]-1)*rows_plhs] = hist_g[0][zz]/norm_rgb;
//                 labeled_element_descriptors[zz+2*n_bins_rgb+((int)labeled_elements[ii]-1)*rows_plhs] = hist_b[0][zz]/norm_rgb;
                

                
            }
            
            for(int zz=0; zz<n_bins_grad; zz++){
                
                labeled_element_descriptors[zz+3*n_bins_rgb+((int)labeled_elements[ii]-1)*rows_plhs] = hist_grad[zz]/norm_grad;
                
//                 labeled_element_descriptors[zz+3*n_bins_rgb+((int)labeled_elements[ii]-1)*rows_plhs] = hist_grad[0][zz]/norm_grad;
                
            }
        }
    }
    
    
    mxDestroyArray(mat_i_R_pr);
    mxDestroyArray(mat_i_G_pr);
    mxDestroyArray(mat_i_B_pr);
    mxDestroyArray(b_rgb_pr);
    mxDestroyArray(b_grad_pr);
    mxDestroyArray(r_surr_pr);
    mxDestroyArray(g_surr_pr);
    mxDestroyArray(b_surr_pr);
    mxDestroyArray(hist_r_pr);
    mxDestroyArray(hist_g_pr);
    mxDestroyArray(hist_b_pr);
    mxDestroyArray(hist_grad_pr);
    
}


    

double *hist_c(double *v, int size_v, double *bin_edges, int n_bins){
 
//     size_t size_v          = v.shape()[0];
//     size_t size_bin_edges  = bin_edges.shape()[0];
//     
    /* Output allocation */    
    
    mxArray *n_pr = mxCreateDoubleMatrix(1,n_bins-1,mxREAL);
    double *n = mxGetPr(n_pr);
    
    for(int i=0; i<size_v; i++)
    {
        for(int j=0; j<n_bins; j++)
        {
            
            if(bin_edges[j-1]<=v[i] && v[i]<bin_edges[j])
                n[j-1] = n[j-1] + 1;

        }
        
        if(v[i] == bin_edges[n_bins-1])
            n[n_bins-2] = n[n_bins-2] + 1;  
    }
    
//     cout << n[0] << endl;
    
    return n;
    
}


int norm(double *hist_array, int n_bins){
    
    int n = 0;
    for(int i=0; i<n_bins; i++){
        
        n = n + hist_array[i];
        
    }
 
    return n;
}

// int norm(MatlabMultiArray<double> hist_array, int n_bins){
//     
//     int n = 0;
//     for(int i=0; i<n_bins; i++){
//         
//         n = n + hist_array[0][i];
//         
//     }
//  
//     return n;
// }



// MatlabMultiArray<double> hist_c(double *v, int size_v, double *bin_edges, int n_bins){
//  
// //     size_t size_v          = v.shape()[0];
// //     size_t size_bin_edges  = bin_edges.shape()[0];
// //     
//     /* Output allocation */    
//     
//     mxArray *n_pr = mxCreateDoubleMatrix(1,n_bins-1,mxREAL);
//     MatlabMultiArray<double> n(n_pr);
//     
//     for(int i=0; i<size_v; i++)
//     {
//         for(int j=0; j<n_bins; j++)
//         {
//             
//             if(bin_edges[j-1]<=v[i] && v[i]<bin_edges[j])
//                 n[0][j-1] = n[0][j-1] + 1;
// 
//         }
//         
//         if(v[i] == bin_edges[n_bins-1])
//             n[0][n_bins-2] = n[0][n_bins-2] + 1;  
//     }
//     
// //     cout << n[0] << endl;
//     
//     return n;
//     
// }



 