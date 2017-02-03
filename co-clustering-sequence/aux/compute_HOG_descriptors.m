function [ descriptors ] = compute_HOG_descriptors( I, labeled_elements )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[Gmag,Gdir] = imgradient(rgb2gray(I));

Gdir = Gdir.*double(Gdir~=(-180)) + 180*double(Gdir==(-180));

step=22.5;
total=0;
count=0;
Gdir_q = 0;

while(total<180)  
    Gdir_q = Gdir_q + (count+1)*double(Gdir>(count*step) & Gdir<=((count+1)*step));
    total = total + step;
    count = count + 1;
end

count2 = count;
count  = 0;
total = -180;

while(total<0)  
    Gdir_q = Gdir_q + (count2+1)*double(Gdir>(-180+(count*step)) & Gdir<=(-180+(count+1)*step));
    total = total + step;
    count = count + 1;
    count2 = count2 + 1;
end


[s_x_i, s_y_i, ~] = size(I);
[s_x, s_y] = size(labeled_elements);
l_window   = 3;
bins       = 16;
pos_bins   = [1:16];

n_elems = sum(sum(labeled_elements~=0));
descriptors = zeros(4*bins, n_elems);

for i=1:s_x
    for j=1:s_y
        if(labeled_elements(i,j)~=0)
            if(mod(i,2)==0)
                x = i/2;
            else
                x = (i-1)/2;
            end
            
            if(mod(j,2)==0)
                y = j/2;
            else
                y = (j-1)/2;
            end
            
                v1    = zeros(1,bins);
                v2    = zeros(1,bins);
                v3    = zeros(1,bins);
                v4    = zeros(1,bins);

                for ii=1:l_window
                    for jj=1:l_window
                        if(x-l_window-1+ii>=1 && y-l_window-1+jj>=1 && x-l_window-1+ii<=s_x_i && y-l_window-1+jj<=s_y_i)
                            v1(Gdir_q(x-l_window-1+ii,y-l_window-1+jj)) = v1(Gdir_q(x-l_window-1+ii,y-l_window-1+jj)) + Gmag(x-l_window-1+ii,y-l_window-1+jj);
                        end
                        
                        if(x-l_window-1+ii>=1 && y-1+jj>=1 && x-l_window-1+ii<=s_x_i && y-1+jj<=s_y_i)
                            v2(Gdir_q(x-l_window-1+ii,y-1+jj))          = v2(Gdir_q(x-l_window-1+ii,y-1+jj)) + Gmag(x-l_window-1+ii,y-1+jj);
                        end
                        
                        if(x-1+ii>=1 && y-1+jj>=1 && x-1+ii<=s_x_i && y-1+jj<=s_y_i)
                            v3(Gdir_q(x-1+ii,y-1+jj))                   = v3(Gdir_q(x-1+ii,y-1+jj))      + Gmag(x-1+ii,y-1+jj);
                        end
                        
                        if(x-1+ii>=1 && y-l_window-1+jj>=1 && x-1+ii<=s_x_i && y-l_window-1+jj<=s_y_i)
                            v4(Gdir_q(x-1+ii,y-l_window-1+jj))          = v4(Gdir_q(x-1+ii,y-l_window-1+jj)) + Gmag(x-1+ii,y-l_window-1+jj);
                        end

                    end
                end    

                if(sum(v1)>0)
                    total_v1 = sum(v1);
                else
                    total_v1 = 1;
                end
                
                if(sum(v2)>0)
                    total_v2 = sum(v2);
                else
                    total_v2 = 1;
                end
                
                if(sum(v3)>0)
                    total_v3 = sum(v3);
                else
                    total_v3 = 1;
                end
                
                if(sum(v4)>0)
                    total_v4 = sum(v4);
                else
                    total_v4 = 1;
                end
                
                descriptors(:, labeled_elements(i,j)) = [v1'/total_v1; v2'/total_v2; v3'/total_v3; v4'/total_v4];
            
        end        
    end
end

end

