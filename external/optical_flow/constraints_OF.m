function [A, b, Aeq, beq] = constraints_OF(regions_i, regions_j, RAG, variables)

%%%% Count number of 3-cliques %%%%
total_regions = regions_i + regions_j;

conditions = 0;
for i=1:total_regions
    for j=(i+1):total_regions
        if(RAG(i,j)==1)
            %dij
            for k=1:total_regions
                if(RAG(j,k)==1 && k~=i)
                    if(RAG(k,i)==1)
                        conditions = conditions +1;
                    end
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = sparse(conditions, variables);
b = zeros(1,conditions);
row = 1;

for i=1:total_regions
    for j=(i+1):total_regions
        if(RAG(i,j)==1)
            %dij
            for k=1:total_regions
                if(RAG(j,k)==1 && k~=i)
                    if(RAG(k,i)==1)
                        if(j>i)
                            A(row,pos_d(i,j,total_regions)) = 1;
                        else
                            A(row,pos_d(j,i,total_regions)) = 1;
                        end
                        if(k>j)
                            A(row,pos_d(j,k,total_regions)) = -1;
                        else
                            A(row,pos_d(k,j,total_regions)) = -1;
                        end
                        if(i>k)
                            A(row,pos_d(k,i,total_regions)) = -1;
                        else
                            A(row,pos_d(i,k,total_regions)) = -1;
                        end
                        
                        row = row + 1;
                    end
                end
            end
        end
    end
end


pos = 0;
for i=1:regions_i
    for j=i+1:regions_i
        pos = pos+1;
    end
end

pos2 = 0;
for i=regions_i+1:regions_i+regions_j
    for j=1:regions_i
    end
    pos2 = pos2 +1;
end

Aeq = zeros(pos+pos2,variables);
beq = ones(1,pos+pos2);

row = 1;
for i=1:regions_i
    for j=i+1:regions_i
        Aeq(row,pos_d(i,j,total_regions)) = 1;
        row = row + 1;
    end
end

for i=regions_i+1:regions_i+regions_j
    for j=1:regions_i
        Aeq(row,pos_d(j,i,total_regions)) = -1;
        beq(row) = -regions_i+1;
    end
    row = row + 1;
end


end