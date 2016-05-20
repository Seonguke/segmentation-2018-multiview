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
#include <map>

using namespace std;    // For std::vector

typedef struct{
    
    double real;
    double imag;
    
}TNum;

complex<double> mycomplex (0.0,1.0);
double PI_trans = 3.1416/180;

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[] )
{
    if(nrhs != 3)
        mexErrMsgTxt("There should be 3 inputs (partition, angles and labeled_elements)");
    
    
    
    double  *partition, *angles, *labeled_elements, *gPb_angles;
    int m_B, n_B, z_B, m_P, n_P, z_P, reg, els, window;
    complex<double> complx_B_scalar;
    
    /* Input parammeters
     * prhs[0] -> partition         (double(matrix))
     * prhs[1] -> angles            (double(matrix))
     * prhs[2] -> labeled_elements  (double(matrix))
     */
    
    partition = mxGetPr(prhs[0]);
    angles = mxGetPr(prhs[1]);
    labeled_elements = mxGetPr(prhs[2]);
    
    m_P = mxGetM(prhs[0]);
    n_P = mxGetN(prhs[0]);
    z_P = m_P * n_P;
    
    m_B = mxGetM(prhs[1]);
    n_B = mxGetN(prhs[1]);
    z_B = m_B * n_B;
    
//     cout << z_B << endl;

    
    mxArray *mat_p_pr_mat = mxCreateDoubleMatrix(m_B,n_B,mxREAL);
    MatlabMultiArray<double> mat_p_mat(mat_p_pr_mat);
    //double *mat_P = mxGetPr(mat_P_pr);
   
    
    int idx_partition = 0;
    
    for(int n=1; n < n_B; n=n+2){
        for(int m=1; m < m_B; m=m+2){
            
            mat_p_mat[m][n] = partition[idx_partition];
            idx_partition++;
            
        }
    }
    
//     plhs[0] = mxCreateDoubleMatrix(m_B,n_B,mxREAL);
//     double *mat_p = mxGetPr(plhs[0]);
    
    mxArray *mat_p_pr = mxCreateDoubleMatrix(m_B,n_B,mxREAL);
    double *mat_p = mxGetPr(mat_p_pr);
    int idx_mat_p = 0;
    
    for(int n=0; n < n_B; n++){
        for(int m=0; m < m_B; m++){
            
            mat_p[idx_mat_p] = mat_p_mat[m][n];
            idx_mat_p++;
            
        }
    }
    
    // Look for number of elements and regions, we will use this iteration to calculate gPb_angles and to
    // tranform the partition to the new space B too
    els = 0;
    reg = 0;
    
    
    plhs[1] = mxCreateDoubleMatrix(m_B,n_B,mxREAL);
    gPb_angles = mxGetPr(plhs[1]);
    
    
    
    for(int i=0; i < z_B; i++){
        
        if(els < labeled_elements[i])
            els = labeled_elements[i];
        
        gPb_angles[i] = (angles[i]-1)*22.5 + (angles[i]==1)*180;
        
    }
    
    for(int i=0; i < z_P; i++){
        
        if(reg < partition[i])
            reg = partition[i];
        
    }
    
    
//      cout << els << " " << reg << endl;
    
    /* Output parammeters
    * plhs[0] -> B              (sparse(matrix))
    * plhs[1] -> gPb_angles     (double(matrix)) 
    */
    
//     mxArray *B_pr = mxCreateDoubleMatrix(els,reg,mxCOMPLEX);
//     double *B_r = mxGetPr(B_pr);
//     double *B_i = mxGetPi(B_pr);
    
//     plhs[0] = mxCreateDoubleMatrix(els,reg,mxCOMPLEX);
//     double *B_r = mxGetPr(plhs[0]);
//     double *B_i = mxGetPi(plhs[0]);
 
    int n_diff_zero = 0;
    vector<int> index;
    vector<TNum> B_data;
    
    
    //MatlabMultiArray complex<double> B(plhs[0]);
    
    window = 5;
    
    for(int ii = 0; ii < z_B; ii++){
     
         if(labeled_elements[ii] != 0){
                 
             if((((ii+1)%m_B)%2) != 0){        // mod(mod(ii,m_l),2)    In c++ is == (from 0 to end-1) in matlab is != (from 1 to end)
             //Horizontal contours
                 
                 //Problem with 180 degrees, search for vertical contours
                 if(gPb_angles[ii] == 180){
                     
                     for(int jj = -window; jj <= window; jj=jj+2){
                         
                         if(labeled_elements[ii-m_B-jj] != 0){   // If a vertical contour is found, compare if the regions to look at are the same
                             
                             if((mat_p[ii-m_B-m_B-jj] != mat_p[ii+m_B-m_B-jj])){
                                 
                                 if(((mat_p[ii-m_B-m_B-jj] == mat_p[ii-1])||(mat_p[ii-m_B-m_B-jj] == mat_p[ii+1]))&&
                                         (mat_p[ii+m_B-m_B-jj] == mat_p[ii-1])||(mat_p[ii+m_B-m_B-jj] == mat_p[ii+1])){
                                     
                                     n_diff_zero ++;
                                     
                                     complx_B_scalar = -exp(mycomplex*(gPb_angles[ii])*PI_trans);
//                                      B_r[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-m_B-m_B-jj]-1)] = complx_B_scalar.real();
//                                      B_i[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-m_B-m_B-jj]-1)] = complx_B_scalar.imag();
                                     index.push_back(((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-m_B-m_B-jj]-1));
                                     B_data.push_back({complx_B_scalar.real(),complx_B_scalar.imag()});
                                     
                                     complx_B_scalar = exp(mycomplex*(gPb_angles[ii])*PI_trans);
//                                      B_r[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B-m_B-jj]-1)] = complx_B_scalar.real();
//                                      B_i[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B-m_B-jj]-1)] = complx_B_scalar.imag();
                                     index.push_back(((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B-m_B-jj]-1));
                                     B_data.push_back({complx_B_scalar.real(),complx_B_scalar.imag()});
                                     
                                     break;
                                 }
                                 
                             }
                         }
                         
                         else if(labeled_elements[ii+m_B-jj] != 0){
                             
                             if((mat_p[ii+m_B-m_B-jj] != mat_p[ii+m_B+m_B-jj])){
                                 
                                 if(((mat_p[ii+m_B-m_B-jj] == mat_p[ii-1])||(mat_p[ii+m_B-m_B-jj] == mat_p[ii+1]))&&
                                         (mat_p[ii+m_B+m_B-jj] == mat_p[ii-1])||(mat_p[ii+m_B+m_B-jj] == mat_p[ii+1])){
                                     
                                     n_diff_zero ++;
                                     
                                     complx_B_scalar = exp(mycomplex*(gPb_angles[ii])*PI_trans);
//                                      B_r[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B+m_B-jj]-1)] = complx_B_scalar.real();
//                                      B_i[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B+m_B-jj]-1)] = complx_B_scalar.imag();
                                     index.push_back(((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B+m_B-jj]-1));
                                     B_data.push_back({complx_B_scalar.real(),complx_B_scalar.imag()});
                                     
                                     
                                     complx_B_scalar = -exp(mycomplex*(gPb_angles[ii])*PI_trans);
//                                      B_r[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B-m_B-jj]-1)] = complx_B_scalar.real();
//                                      B_i[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B-m_B-jj]-1)] = complx_B_scalar.imag();
                                     index.push_back(((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B-m_B-jj]-1));
                                     B_data.push_back({complx_B_scalar.real(),complx_B_scalar.imag()});
                                     
                                     break;
                                     
                                 }
                                 
                             }
                             
                         }
                         
                         if(jj == window){   // The program hasn't found any matches
                             
                             complx_B_scalar = -exp(mycomplex*(gPb_angles[ii])*PI_trans);
//                               B_r[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-1]-1)] = complx_B_scalar.real();
//                               B_i[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-1]-1)] = complx_B_scalar.imag();
                             index.push_back(((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-1]-1));
                             B_data.push_back({complx_B_scalar.real(),complx_B_scalar.imag()});
                             
                             complx_B_scalar = exp(mycomplex*(gPb_angles[ii])*PI_trans);
//                               B_r[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+1]-1)] = complx_B_scalar.real();
//                               B_i[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+1]-1)] = complx_B_scalar.imag();
                             index.push_back(((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+1]-1));
                             B_data.push_back({complx_B_scalar.real(),complx_B_scalar.imag()});
                             
                             n_diff_zero ++;
                         }
                         
                     }
                 }
                 else{
                     
                     complx_B_scalar = -exp(mycomplex*(gPb_angles[ii])*PI_trans);
//                      B_r[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-1]-1)] = complx_B_scalar.real();
//                      B_i[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-1]-1)] = complx_B_scalar.imag();
                     index.push_back(((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-1]-1));
                     B_data.push_back({complx_B_scalar.real(),complx_B_scalar.imag()});
                     
                     complx_B_scalar = exp(mycomplex*(gPb_angles[ii])*PI_trans);
//                      B_r[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+1]-1)] = complx_B_scalar.real();
//                      B_i[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+1]-1)] = complx_B_scalar.imag();
                     index.push_back(((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+1]-1));
                     B_data.push_back({complx_B_scalar.real(),complx_B_scalar.imag()});
                     
                     n_diff_zero ++;
                 }
             }
             else{
                 //Vertical contours
                 
                 //Problem with 90 degrees, search for vertical contours
                 if(gPb_angles[ii] == 90){
                     
                     for(int jj = -window; jj <= window; jj++){
                         
                         if(labeled_elements[ii-m_B-jj] != 0){   // If an horizontal contour is found, compare if the regions to look at are the same
                             
                             if((mat_p[ii-m_B-1-jj] != mat_p[ii-m_B+1-jj])){
                                 if(((mat_p[ii-m_B-1-jj] == mat_p[ii-1])||(mat_p[ii-m_B+1-jj] == mat_p[ii+1]))&&
                                         (mat_p[ii-m_B+1-jj] == mat_p[ii-1])||(mat_p[ii-m_B-1-jj] == mat_p[ii+1])){
                                     
                                     n_diff_zero ++;
                                     
                                     complx_B_scalar = -exp(mycomplex*(gPb_angles[ii])*PI_trans);
//                                      B_r[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-m_B-1-jj]-1)] = complx_B_scalar.real();
//                                      B_i[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-m_B+1-jj]-1)] = complx_B_scalar.imag();
                                     index.push_back(((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-m_B+1-jj]-1));
                                     B_data.push_back({complx_B_scalar.real(),complx_B_scalar.imag()});
                                     
                                     complx_B_scalar = exp(mycomplex*(gPb_angles[ii])*PI_trans);
//                                      B_r[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-m_B-1-jj]-1)] = complx_B_scalar.real();
//                                      B_i[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-m_B+1-jj]-1)] = complx_B_scalar.imag();
                                     index.push_back(((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-m_B+1-jj]-1));
                                     B_data.push_back({complx_B_scalar.real(),complx_B_scalar.imag()});
                                     
                                     break;
                                 }
                             }
                         }
                         else if(labeled_elements[ii+m_B-jj] != 0){
                             
                             if((mat_p[ii+m_B-1-jj] != mat_p[ii+m_B+1-jj])){
                                 if(((mat_p[ii+m_B-1-jj] == mat_p[ii-1])||(mat_p[ii+m_B+1-jj] == mat_p[ii+1]))&&
                                         (mat_p[ii+m_B+1-jj] == mat_p[ii-1])||(mat_p[ii+m_B-1-jj] == mat_p[ii+1])){
                                     
                                     n_diff_zero ++;
                                     
                                     complx_B_scalar = exp(mycomplex*(gPb_angles[ii])*PI_trans);
//                                      B_r[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B-1-jj]-1)] = complx_B_scalar.real();
//                                      B_i[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B+1-jj]-1)] = complx_B_scalar.imag();
                                     index.push_back(((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B+1-jj]-1));
                                     B_data.push_back({complx_B_scalar.real(),complx_B_scalar.imag()});
                                     
                                     complx_B_scalar = -exp(mycomplex*(gPb_angles[ii])*PI_trans);
//                                      B_r[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B-1-jj]-1)] = complx_B_scalar.real();
//                                      B_i[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B+1-jj]-1)] = complx_B_scalar.imag();
                                     index.push_back(((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B+1-jj]-1));
                                     B_data.push_back({complx_B_scalar.real(),complx_B_scalar.imag()});
                                     
                                     break;
                                     
                                 }
                             }
                             
                         }
                     }
                 }
                 else if(gPb_angles[ii] < 90){
                     
                     complx_B_scalar = exp(mycomplex*(gPb_angles[ii])*PI_trans);
//                      B_r[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-m_B]-1)] = complx_B_scalar.real();
//                      B_i[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-m_B]-1)] = complx_B_scalar.imag();
                     index.push_back(((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-m_B]-1));
                     B_data.push_back({complx_B_scalar.real(),complx_B_scalar.imag()});
                     
                     complx_B_scalar = -exp(mycomplex*(gPb_angles[ii])*PI_trans);
//                      B_r[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B]-1)] = complx_B_scalar.real();
//                      B_i[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B]-1)] = complx_B_scalar.imag();
                     index.push_back(((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B]-1));
                     B_data.push_back({complx_B_scalar.real(),complx_B_scalar.imag()});
                     
                     n_diff_zero ++;
                 }
                 else{
                     
                     
                     
                     complx_B_scalar = -exp(mycomplex*(gPb_angles[ii])*PI_trans);
//                      B_r[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-m_B]-1)] = complx_B_scalar.real();
//                      B_i[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-m_B]-1)] = complx_B_scalar.imag();
                     index.push_back(((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii-m_B]-1));
                     B_data.push_back({complx_B_scalar.real(),complx_B_scalar.imag()});
                     
                     complx_B_scalar = exp(mycomplex*(gPb_angles[ii])*PI_trans);
//                      B_r[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B]-1)] = complx_B_scalar.real();
//                      B_i[((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B]-1)] = complx_B_scalar.imag();
                     index.push_back(((int)labeled_elements[ii]-1)+els * ((int)mat_p[ii+m_B]-1));
                     B_data.push_back({complx_B_scalar.real(),complx_B_scalar.imag()});
                     
                     n_diff_zero ++;
                 }
             }
         }
    }

    
//     sort(index.begin(),index.end());
    
    multimap<int,TNum> map;

    for(int i = 0; i < index.size(); i++){
        
        map.insert(pair<int,TNum> (index[i],B_data[i]));    
        
    }
   
    int idx_map = 0;
    
    for (multimap<int, TNum>::iterator it = map.begin(); it != map.end(); it++){
        
        index[idx_map] = (*it).first;
        B_data[idx_map].real = (*it).second.real;
        B_data[idx_map].imag = (*it).second.imag;
        
        idx_map++;
                
//      cout << "  [" << (*it).first << ", " << (*it).second.real << "+" << (*it).second.imag << "]" << endl;
        
    }
    
    
//     cout << n_diff_zero*2 << endl;
    
//     plhs[0] = mxCreateDoubleMatrix(index.size(),1,mxREAL);
//     double *out = mxGetPr(plhs[0]);
//
//     for(int i = 0; i < index.size(); i++)
//         out[i] = index[i]+1;
    
    plhs[0] = mxCreateSparse(els,reg,n_diff_zero*2,mxCOMPLEX);
    
    mwIndex *irs, *jcs;
    double *sr, *si;
    
    sr = mxGetPr(plhs[0]);
    si = mxGetPi(plhs[0]);
    irs = mxGetIr(plhs[0]);
    jcs = mxGetJc(plhs[0]);
    
    int idx, accB;
    idx = 0;
    accB = 0;
    
    float r;
    r = 0;
    
    jcs[idx] = 0;
    idx ++;
    
//      cout << index.size() << endl;
//      cout << n_diff_zero*2 << endl;
    
    for(int i = 0; i < index.size(); i++){
        
        if(i != 0){
            
            if(((index[i]/els)%reg + 1) > ((index[i-1]/els)%reg + 1)){
                
                jcs[idx] = accB;
                idx ++;
            }
            
            if(i == index.size()-1){
                
                jcs[idx] = accB+1;
                idx ++;
            }
            
            // Case with a regions without labeled elements (example: 1 pixel region)
            
            if(i != index.size()-1){
                
                if((((index[i]/els)%reg + 1)-((index[i-1]/els)%reg + 1)) > 1){
                    
                    r = (index[i] - index[i-1])/els;        // If two or more regions without labels are consecutive
                    r = floor(r);
                    
                    for(int pass = 0; pass < r; pass ++){
                        
                        jcs[idx] = jcs[idx-1];
                        idx++;
                        
                    }
                }
            }
        }
        
        // Case with first region without labels
        
        else if((i==0) && (index[i] > els)){        
            
            r = index[i]/els;        // If two or more regions without labels are consecutive
            r = floor(r);
            
            for(int pass = 0; pass < r; pass ++){
                
                jcs[idx] = accB;
                idx++;
                
            }
        }
        
        sr[accB] = B_data[i].real;
        si[accB] = B_data[i].imag;
        irs[accB] = index[i]%els;
        accB ++;
        
        
    }

    mxDestroyArray(mat_p_pr_mat);
    mxDestroyArray(mat_p_pr);
//     mxDestroyArray(B_pr);
    
}