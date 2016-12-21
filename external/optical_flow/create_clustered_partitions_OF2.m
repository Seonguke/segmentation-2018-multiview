function [partition_i, partition_j] = create_clustered_partitions_OF2(I_i, I_j, seg_i, seg_j, RAG, x)
regions_i = max(max(seg_i));
regions_j = max(max(seg_j));

total_regions = regions_i + regions_j;


partition = x2partition(x', total_regions);

R1 = I_i(:,:,1);
G1 = I_i(:,:,2);
B1 = I_i(:,:,3);

R2 = I_j(:,:,1);
G2 = I_j(:,:,2);
B2 = I_j(:,:,3);

color1 = zeros(regions_i,3);
for i=1:regions_i
    color1(i,:) = [mean(R1(seg_i==i)) mean(G1(seg_i==i)) mean(B1(seg_i==i))];
end

color2 = zeros(regions_j,3);
for i=1:regions_j
    color2(i,:) = [mean(R2(seg_j==i)) mean(G2(seg_j==i)) mean(B2(seg_j==i))];
end

color = [color1; color2];
[sx sy ch] = size(I_i);

tmp = jet(total_regions);
color = uint8(255*tmp(randperm(total_regions,total_regions),:));

partition1 = zeros(sx, sy, ch);
partition2 = zeros(sx, sy, ch);

color3 = zeros(1,total_regions);
c_index=1;

partition = partition.*(RAG + eye(total_regions));

%%%% Label regions from seg_i %%%%
for i=1:regions_i
    color3(i) = i;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Label regions from seg_j %%%%
%%%%   with labels of seg_i   %%%%
for i=(regions_i+1):total_regions
    assigned_labels = partition(i,1:regions_i);
    %%%% Analyze labels assigned to region %%%%
    if(sum(assigned_labels) == 1)
        %% Label directly assigned from seg_i %%
        color3(i) = find(assigned_labels==1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=(regions_i+1):total_regions
    if(color3(i)==0)
        %%  If no label is assigned from seg_i, %%
        %%  assign most likely neightbour label %%
            clusters = partition(i,:);
            regions = find(clusters);
            
            label_neightbors = zeros(1,regions_i); 
            for j=1:length(regions)
                if(color3(regions(j))~=0)
                    label_neightbors(color3(regions(j))) = label_neightbors(color3(regions(j))) + 1; 
                end
            end
            
            max_label = 0;
            max_value = -1;
            for k=1:length(regions_i)
                if(label_neightbors(k)>max_value)
                    max_label = k;
                    max_value = label_neightbors(k);
                end
            end
            
            color3(i) = max_label;
            
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:total_regions 
%     if(sum(assigned_labels(i,:))==0) 
%         1 
%     end
% end


% for i=(regions_i+1):total_regions
%     clusters = partition(i,:);
%     regions = find(clusters);
%     label = 0;
% 
%     for j=1:regions_i
%         if(clusters(j) == 1 & RAG(i,j)==1)
%             label = j;
%         end
%     end
%     
%     if(i==328)
%         1
%     end
%     if(label==0)
%     %%% Search labels from neightbours %%%
%     for j=1:length(regions)
%          if(regions(j)==331)
%                 2
%             end
%         if(color3(regions(j))~=0 & RAG(i,regions(j))==1)
%             label = color3(regions(j));
%         end
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     end
%     color3(i) = label;
% 
% end

% for i=1:regions_i
%     clusters = partition(i,:);
%     regions = find(clusters);
%     label = 0;
% 
%     %%% Search labels from neightbours %%%
%     for j=1:length(regions)
%         if(color3(regions(j))~=0 & RAG(i,regions(j))==1)
%             label = color3(regions(j));
%         end
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     if (label==0)
%         label = c_index;
%     end
% 
%     for j=1:length(regions)
%         if(color3(regions(j))==0 & RAG(i,regions(j))==1)
%             color3(regions(j)) = label;
%         end
%     end
% 
%     color3(i) = label;
%     c_index=c_index+1;
% end


% for i=regions_i+1:total_regions
%     if(color3(i)==0)
%         color3(i) = c_index;
%         c_index=c_index+1;
%     end
% end

s1 = size(seg_i,1);
s2 = size(seg_i,2);

partition_i = zeros(s1,s2);
partition_j = zeros(s1,s2);

for i=1:s1
    for j=1:s2
        partition_i(i,j) = color3(seg_i(i,j));
        partition_j(i,j) = color3(regions_i+seg_j(i,j));
    end
end

% for i=1:total_regions
% %     clusters = partition(i,:);
% %     region = 0;
% %     for j=1:total_regions
% %         if(clusters(j)==1 && j~=i)
% %             if(region == 0)
% %                 region = j;
% %             else
% %                 region = -1;
% %             end
% %         end
% %     end
% %     if(region==0)
% %         region = i;
% %     end
%     
%     if(i<=regions_I1)
%         mask = uint8(seg1==i);
% %         if(region > 0)
%             mask_r = mask*color(color3(i),1);
%             mask_g = mask*color(color3(i),2);
%             mask_b = mask*color(color3(i),3);
% %         else
% %             mask_r = mask*255;
% %             mask_g = mask*0;
% %             mask_b = mask*0;
% %         end 
%         partition1(:,:,1) = uint8(partition1(:,:,1)) + mask_r;
%         partition1(:,:,2) = uint8(partition1(:,:,2)) + mask_g;
%         partition1(:,:,3) = uint8(partition1(:,:,3)) + mask_b;
%     else
%         mask = uint8(seg2==(i-regions_I1));
% %         if(region > 0)
%             mask_r = mask*color(color3(i),1);
%             mask_g = mask*color(color3(i),2);
%             mask_b = mask*color(color3(i),3);
% %         else
% %             mask_r = mask*255;
% %             mask_g = mask*0;
% %             mask_b = mask*0;
% %         end 
%         partition2(:,:,1) = uint8(partition2(:,:,1)) + mask_r;
%         partition2(:,:,2) = uint8(partition2(:,:,2)) + mask_g;
%         partition2(:,:,3) = uint8(partition2(:,:,3)) + mask_b;
%     end
%     
% end
% 
% figure; imshow(uint8(partition1));
% figure; imshow(uint8(partition2));
% 
% % figure; imshow(I1)
end