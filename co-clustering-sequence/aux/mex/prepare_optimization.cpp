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
    ConstMatlabMultiArray<double> regions_adjacency_total(prhs[0]);
    ConstMatlabMultiArray<double> partitions_adjacency_total(prhs[1]);

    variables = mxGetN(prhs[0]);
    
    ConstMatlabMultiArray<double> regions_partitions(prhs[2]);
    int num_partitions = mxGetN(prhs[2]);
    
    double total_regions = 0;
    MatlabMultiArray<double> cumulate = mxCreateDoubleMatrix(1, num_partitions+1, mxREAL);

    for(int ii=0; ii<num_partitions; ii++)
    {
        total_regions = total_regions + regions_partitions[0][ii];
        cumulate[0][ii+1] = total_regions;
    }
    
    std::cout << "Variables: " << variables << std::endl;
    std::cout << "Num partitions: " << num_partitions << std::endl;
    std::cout << "Total regions: " << total_regions << std::endl;
    std::cout << "Cumulate: " << cumulate[0][2] << std::endl;
        
    /* Output allocation */    
    plhs[0] = mxCreateDoubleMatrix(2, variables, mxREAL);
    MatlabMultiArray<double> regions_adjacency_total_offset(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(total_regions, total_regions, mxREAL);
    MatlabMultiArray<double> regions_variables_offset(plhs[1]);
    
    plhs[2] = mxCreateDoubleMatrix(1, num_partitions+1, mxREAL);
    MatlabMultiArray<double> cumulate_out(plhs[2]);
      
    cumulate_out = cumulate;
    
    for(int ii=0; ii<variables; ii++)
    {
        regions_adjacency_total_offset[0][ii] = regions_adjacency_total[0][ii] + cumulate[0][partitions_adjacency_total[0][ii]-1];
        regions_adjacency_total_offset[1][ii] = regions_adjacency_total[1][ii] + cumulate[0][partitions_adjacency_total[1][ii]-1];

        regions_variables_offset[regions_adjacency_total_offset[0][ii]-1][regions_adjacency_total_offset[1][ii]-1] = ii+1;
        regions_variables_offset[regions_adjacency_total_offset[1][ii]-1][regions_adjacency_total_offset[0][ii]-1] = ii+1;
    }
    
}