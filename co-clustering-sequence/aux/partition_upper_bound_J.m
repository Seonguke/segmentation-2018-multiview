function [ J, partition ] = partition_upper_bound_J( seg, ground_truth, N_regions, fixed_regions)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if(nargin==2)
    N_regions     = 0;
    fixed_regions = 0;
elseif(nargin==3)
    fixed_regions = 0;
end

%% Variables allocation
%regions_i = 2;
regions_seg = max(max(seg));

[s_x, s_y] = size(ground_truth);
intersection = zeros(regions_seg, 1);
delta_union = zeros(regions_seg, 1);
area_mask = sum(sum(ground_truth~=0));

%% Cost and weight functions
for i=1:s_x
    for j=1:s_y
        if(ground_truth(i,j)~=0)
            intersection(seg(i,j)) = intersection(seg(i,j)) + 1;
        else
            delta_union(seg(i,j)) = delta_union(seg(i,j)) + 1;
        end
    end
end

t = [intersection' 0];
f = [delta_union' area_mask];

%% Equality restrictions
Aeq = [ zeros(1,length(intersection)) 1];
Beq = 1;

%% Equality restrictions - number of regions fixed
%N_regions = 0;
if(N_regions~=0)
    if(fixed_regions == 1)
        Aeq = [Aeq;
           ones(1,length(intersection)) 0];
        Beq = [Beq;
           N_regions];
       
       A = zeros(1,length(intersection)+1);
       B = 0;
       
    else
        A = [ones(1,length(intersection)) 0];
        B = N_regions;
    end
else
    A = zeros(1,length(intersection)+1);
    B = 0;
end

% %% Equality restrictions - number of regions fixed
% %N_regions = 0;
% if(N_regions~=0)
%     Aeq = [Aeq;
%            ones(1,length(intersection)) 0];
%     Beq = [Beq;
%            N_regions];
% end
% 
% %% Inequality restrictions 
% A = zeros(1,length(intersection)+1);
% B = 0;

% \epsilon value considered zero and maximum number of iterations
threshold = 1e-5;
num_iter_max = 1000;

% Initial values
delta_bar = 0;
h_delta_bar = threshold + 1;
num_iter = 0;

opt = cplexoptimset('cplex');
    
%% Algorithm to solve an LFCO
while(h_delta_bar>threshold && num_iter<num_iter_max)

    % Linear optimization problem (minus to maximize)
    func = -(t-delta_bar*f);
    [x, fval] = cplexbilp(func',A,B,Aeq,Beq,[], opt);

    h_delta_bar = -fval;

    % Update delta_bar
    delta_bar = (t*x)/(f*x);
    num_iter = num_iter + 1;
end
    
[ partition ] = visualize_jaccard_partition( x, seg );
J = (t*x)/(f*x);
% figure;imshow(label2rgb(partition))

end

