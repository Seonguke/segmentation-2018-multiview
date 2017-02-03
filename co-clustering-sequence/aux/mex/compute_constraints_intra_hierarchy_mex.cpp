#include "mex.h"
#include "matlab_multiarray.hpp"
#include <iostream>
#include <list>
#include <set>
#include <algorithm>

void mexFunction( int nlhs, mxArray *plhs[], 
        		  int nrhs, const mxArray*prhs[] )
{

    double offset;
    int variables_opt;
    
    if(nrhs<1)
        mexErrMsgTxt("There should be at least 1 input parameter");
    
    /* Input parameters */
    ConstMatlabMultiArray<double> partition(prhs[0]);
    ConstMatlabMultiArray<double> merging_sequence(prhs[1]);
    ConstMatlabMultiArray<double> nodes_leaves(prhs[2]);
    ConstMatlabMultiArray<double> RAG(prhs[3]);
    ConstMatlabMultiArray<double> regions_variables(prhs[4]);
    variables_opt = mxGetScalar(prhs[5]);
    offset = mxGetScalar(prhs[6]);
    
    std::size_t n_fusions = merging_sequence.shape()[0];
    
    std::size_t s_x = partition.shape()[0];
    std::size_t s_y = partition.shape()[1];
    
    int n_regions = 0;
    for(int i=0; i<s_x; i++)
    {
        for(int j=0; j<s_y; j++)
        {
            if(partition[i][j] > n_regions)
                n_regions = partition[i][j];
        }
    }
        
    int last_inter;
    double reg_A, reg_B, reg_C, count_inter, count_intra;
    double leaves_A[n_regions], leaves_B[n_regions];
    double n_total_equalities = 0;
    double n_total_inequalities = 0;
    
    for(int i=0; i < n_fusions-1; i++)
    {
        reg_A = merging_sequence[i][0];
        reg_B = merging_sequence[i][1];
        reg_C = merging_sequence[i][2];
        
        if(reg_A > n_regions | reg_B > n_regions)
        {
            for(int j=0; j < n_regions; j++)
            {
                leaves_A[j] = nodes_leaves[j][reg_A-1];
                leaves_B[j] = nodes_leaves[j][reg_B-1];
            }
            
            last_inter = 0;
            
            count_inter = 0;
            count_intra = 0;
            
            for(int j=0; j < n_regions; j++)
            {
                if(leaves_A[j] == 1)
                {
                    for(int k=0; k < n_regions; k++)
                    {
                        if(leaves_B[k] == 1)
                        {
                            if(RAG[j][k] == 1)
                            {
                                count_inter++;
                            }   
                        }
                        
                        if(leaves_A[k] == 1 & k>j)
                        {
                            if(RAG[j][k] == 1)
                            {
                                count_intra++;
                            }   
                        }
                        
                    }
                }
                
                if(leaves_B[j]==1)
                {
                    for(int k=0; k<n_regions; k++)
                    {
                        if(leaves_B[k]==1 && k>j)
                        {
                            if(RAG[j][k]==1)
                            {
                                count_intra++;
                            }
                            
                        }  
                    }      
                }
            }
            
            if(count_inter > 1)
            {
                n_total_equalities++;
            }

            if(count_intra > 0)
            {
                n_total_inequalities++;
            }

        }
        
    }
    
    //This is for the last equality: the two 
    //last nodes are never fused
    n_total_equalities++;
    
//     std::cout << "N total_equalities: " << n_total_equalities << std::endl;
//     std::cout << "N total_inequalities: " << n_total_inequalities << std::endl;

    /* Output allocation */    
    plhs[0] = mxCreateDoubleMatrix(n_total_equalities, variables_opt, mxREAL);
    MatlabMultiArray<double> total_equalities(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(n_total_inequalities, variables_opt, mxREAL);
    MatlabMultiArray<double> total_inequalities(plhs[1]);
    
    plhs[2] = mxCreateDoubleMatrix(1, n_total_equalities, mxREAL);
    MatlabMultiArray<double> b_equalities(plhs[2]);
    
    plhs[3] = mxCreateDoubleMatrix(1, n_total_inequalities, mxREAL);
    MatlabMultiArray<double> b_inequalities(plhs[3]);
    
    double count_equalities = 0;
    double count_inequalities = 0;
    double inter_subnodes[variables_opt];
    double intra_subnodes[variables_opt];
    
    for(int i=0; i < n_fusions-1; i++)
    {
        reg_A = merging_sequence[i][0];
        reg_B = merging_sequence[i][1];
        reg_C = merging_sequence[i][2];
        
        if(reg_A > n_regions | reg_B > n_regions)
        {
            for(int j=0; j < n_regions; j++)
            {
                leaves_A[j] = nodes_leaves[j][reg_A-1];
                leaves_B[j] = nodes_leaves[j][reg_B-1];
            }
            
            for(int j=0; j < variables_opt; j++)
            {
                inter_subnodes[j] = 0;
                intra_subnodes[j] = 0;
            }
            
            last_inter = 0;
            
            count_inter = 0;
            count_intra = 0;
            
            for(int j=0; j < n_regions; j++)
            {
                if(leaves_A[j] == 1)
                {
                    for(int k=0; k < n_regions; k++)
                    {
                        if(leaves_B[k] == 1)
                        {
                            if(RAG[j][k] == 1)
                            {
                                inter_subnodes[int(regions_variables[j+offset][k+offset])-1] = 1;
                                count_inter++;
                                last_inter  = regions_variables[j+offset][k+offset]; 
                            }   
                        }
                        
                        if(leaves_A[k] == 1 & k>j)
                        {
                            if(RAG[j][k] == 1)
                            {
                                intra_subnodes[int(regions_variables[j+offset][k+offset])-1] = 1;
                                count_intra++;
                            }   
                        }
                    }
                }
                
                if(leaves_B[j]==1)
                {
                    for(int k=0; k<n_regions; k++)
                    {
                        if(leaves_B[k]==1 && k>j)
                        {
                            if(RAG[j][k]==1)
                            {
                                intra_subnodes[int(regions_variables[j+offset][k+offset])-1] = 1;
                                count_intra++;
                            }
                        }  
                    }      
                }
            }
            
            if(count_inter > 1)
            {
                inter_subnodes[last_inter-1] = 1-count_inter;
                for(int j=0; j < variables_opt; j++)
                {
                    total_equalities[count_equalities][j] = inter_subnodes[j];
                }
                count_equalities++;
            }

            if(count_intra > 0)
            {
                intra_subnodes[last_inter-1] = -count_intra;
                for(int j=0; j < variables_opt; j++)
                {
                    total_inequalities[count_inequalities][j] = intra_subnodes[j];
                }
                count_inequalities++;
            }
        }
    }
    
    // The last two nodes will never be fused
    reg_A = merging_sequence[n_fusions-1][0];
    reg_B = merging_sequence[n_fusions-1][1];
        
    for(int j=0; j < n_regions; j++)
    {
        leaves_A[j] = nodes_leaves[j][reg_A-1];
        leaves_B[j] = nodes_leaves[j][reg_B-1];
    }
    
    double no_fuse_subnodes[variables_opt];
    int count_no_fuse = 0;
    
    for(int j=0; j < variables_opt; j++)
    {
        no_fuse_subnodes[j] = 0;
    }
    
    for(int j=0; j < n_regions; j++)
    {
        if(leaves_A[j] == 1)
        {
            for(int k=0; k < n_regions; k++)
            {
                if(leaves_B[k] == 1)
                {
                    if(RAG[j][k] == 1)
                    {
                        no_fuse_subnodes[int(regions_variables[j+offset][k+offset])-1] = 1;
                        count_no_fuse++;
                    }
                }
            }
        }
    }

    for(int j=0; j < variables_opt; j++)
    {
        total_equalities[count_equalities][j] = no_fuse_subnodes[j];
    }
    b_equalities[0][n_total_equalities-1] = count_no_fuse;                            
    
    
//     % The last two nodes will never be fused
//            
// b_equalities = [b_equalities count_no_fuse];

    
}