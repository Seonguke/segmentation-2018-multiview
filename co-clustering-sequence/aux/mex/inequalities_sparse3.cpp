/*
 * =============================================================
 * fulltosparse.c
 * This example demonstrates how to populate a sparse
 * matrix.  For the purpose of this example, you must pass in a
 * non-sparse 2-dimensional argument of type double.
 *
 * Comment: You might want to modify this MEX-file so that you can 
 * use it to read large sparse data sets into MATLAB.
 *
 * This is a MEX-file for MATLAB.  
 * Copyright (c) 1984-2000 The MathWorks, Inc.
 * =============================================================
 */

/* $Revision: 1.5 $ */

#include <math.h> /* Needed for the ceil() prototype. */
#include "mex.h"
#include <iostream>
#include "matlab_multiarray.hpp"


/* If you are using a compiler that equates NaN to be zero, you 
 * must compile this example using the flag  -DNAN_EQUALS_ZERO.
 * For example:
 *
 *     mex -DNAN_EQUALS_ZERO fulltosparse.c
 *
 * This will correctly define the IsNonZero macro for your C
 * compiler.
 */

#if defined(NAN_EQUALS_ZERO)
#define IsNonZero(d) ((d) != 0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d) != 0.0)
#endif

void mexFunction(
        int nlhs,       mxArray *plhs[],
        int nrhs, const mxArray *prhs[]
        )
{

  /* Check for proper number of input and output arguments. */    
//   if (nrhs != 5) {
//     mexErrMsgTxt("Five input arguments required.");
//   } 
//   if (nlhs > 5) {
//     mexErrMsgTxt("Too many output arguments.");
//   }

  /* Check data type of input argument. */
//   if (!(mxIsDouble(prhs[0]))) {
//     mexErrMsgTxt("Input argument 1 must be of type double.");
//   } 
//   if (!(mxIsDouble(prhs[1]))) {
//     mexErrMsgTxt("Input argument 2 must be of type double.");
//   }
//   if (!(mxIsDouble(prhs[2]))) {
//     mexErrMsgTxt("Input argument 3 must be of type double.");
//   }

//   if (mxGetNumberOfDimensions(prhs[3]) != 2) {
//     mexErrMsgTxt("Input argument 4 must be two dimensional\n");
//   }

    /* Input parameters */
    double total_regions, variables;
    
    /* Input parameters */
    total_regions = mxGetScalar(prhs[0]);
    variables = mxGetScalar(prhs[1]);
    
    ConstMatlabMultiArray<double> regions_adjacency(prhs[2]);
    ConstMatlabMultiArray<double> regions_variables(prhs[3]);
    ConstMatlabMultiArray<double> partitions_adjacency(prhs[4]);
    
    double region_i, region_j;
    double conditions=0;
    double conditions_debug=0;
    
    for(int i=0; i<variables; i++)
    {
        region_i = regions_adjacency[0][i];
        region_j = regions_adjacency[1][i];

        for(int k=0; k<total_regions; k++)
        {
            if(regions_variables[k][region_j-1]!=0 && (k+1)!=region_i)
            {
                conditions_debug++;
                if(regions_variables[k][region_i-1]!=0)
                {
                    conditions++;
                }
            }
        }
        
    }
      
    int rows_matrix = 3*conditions;
         
    /* Allocate output variables */
    plhs[0] = mxCreateSparse(variables, rows_matrix, 3*rows_matrix, mxREAL);
    
    plhs[1] = mxCreateDoubleMatrix(1, rows_matrix, mxREAL);
    MatlabMultiArray<double> b(plhs[1]);
    
    mwIndex *irs, *jcs;
    double *sr, *si;
    
    sr  = mxGetPr(plhs[0]);
    si  = mxGetPi(plhs[0]);
    irs = mxGetIr(plhs[0]);
    jcs = mxGetJc(plhs[0]);  
    
    int acc = 0; 
    int column = 1;
    
    jcs[0] = acc;
    for(int i=0; i<variables; i++)
    {
        region_i = regions_adjacency[0][i];
        region_j = regions_adjacency[1][i];

        for(int k=0; k<total_regions; k++)
        {
            if(regions_variables[k][region_j-1]!=0 && (k+1)!=region_i)
            {
                if(regions_variables[k][region_i-1]!=0)
                {
                    // First column
                    sr[acc] = 1;
                    irs[acc] = regions_variables[region_i-1][region_j-1]-1;
                    acc++;

                    sr[acc] = -1;
                    irs[acc] = regions_variables[region_i-1][k]-1;
                    acc++;

                    sr[acc] = -1;
                    irs[acc] = regions_variables[k][region_j-1]-1;
                    acc++;

                    jcs[column] = acc;
                    column++;

                    // Second column
                    sr[acc] = -1;
                    irs[acc] = regions_variables[region_i-1][region_j-1]-1;
                    acc++;

                    sr[acc] = 1;
                    irs[acc] = regions_variables[region_i-1][k]-1;
                    acc++;

                    sr[acc] = -1;
                    irs[acc] = regions_variables[region_j-1][k]-1;
                    acc++;

                    jcs[column] = acc;
                    column++;

                    // Third column
                    sr[acc] = -1;
                    irs[acc] = regions_variables[region_i-1][region_j-1]-1;
                    acc++;

                    sr[acc] = 1;
                    irs[acc] = regions_variables[k][region_j-1]-1;
                    acc++;

                    sr[acc] = -1;
                    irs[acc] = regions_variables[k][region_i-1]-1;
                    acc++;

                    jcs[column] = acc;
                    column++;
                }         
            }
        }
        
    }  

    for(int i=column; i<(3*conditions+1); i++)
    {
        jcs[i] =acc;
    }
}