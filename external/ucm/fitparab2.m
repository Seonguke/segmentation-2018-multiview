function a = fitparab2(z,ra,rb,theta,filt)
% function [a,b,c] = fitparab2(z,ra,rb,theta)
% % Note: Changed by jpont from fitparab to fitparab2 to avoid crash with 
% %       the same function in segbench/Gradients
%
% Fit cylindrical parabolas to elliptical patches of z at each
% pixel.  
%
% INPUT
%	z	Values to fit.
%	ra,rb	Radius of elliptical neighborhood, ra=major axis.
%	theta	Orientation of fit (i.e. of minor axis).
%
% OUTPUT
%	a,b,c	Coefficients of fit: a + bx + cx^2
%


% compute the interior quickly with convolutions
a = conv2(z,filt(:,:,1),'same');
%fix border with mex file
a = savgol_border(a, z, ra, rb, theta);

