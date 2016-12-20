function [out, V] = MF_analyze(expr, vols)
% _
% Analyze Volume or Volumes
% FORMAT [out, V] = MF_analyze(expr, vols)
% 
%     expr - a string describing a mathematical operation on V
%     vols - an n x 1 cell array specifying file paths of images
% 
%     out  - the result of the calculation described by expr
%     V    - an n x v matrix with all volumes at all voxels
% 
% FORMAT out = MF_analyze(expr, vols) loads all the images specified by
% vols, performs the operation indicated by expr and returns the result.
% The string in expr must contain the capital letter V. All images need
% to have the same dimensions.
% 
% FORMAT [out, V] = MF_analyze(expr, vols) also returns the image matrix.
% 
% Exemplary usage:
%     f_mean = MF_analyze('mean(V)', {'C:\f0001.img'})
% 
% Similar tools:
%     MF_correlate, MF_histogram
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 27/11/2014, 17:50 (V0.2/V8)
%  Last edit: 17/08/2015, 22:30 (V0.3/V11)


% Get expression if required
%-------------------------------------------------------------------------%
if nargin < 1 || isempty(expr)
    expr = spm_input('Expression:', 1, 's', 'mean(V)');
end;

% Get paths if required
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(vols)
    vols = cellstr(spm_select(Inf, 'image', 'Please select all images that you want to analyze!', [], pwd, '.*', '1'));
end;

% Get matrix dimensions
%-------------------------------------------------------------------------%
v1_hdr = spm_vol(vols{1});
v = prod(v1_hdr.dim);
n = length(vols);
clear v1_hdr

% Load all volumes
%-------------------------------------------------------------------------%
V = zeros(n,v);
for i = 1:n
    v_hdr  = spm_vol(vols{i});
    v_img  = spm_read_vols(v_hdr);
    v_img  = reshape(v_img,[1 v]);
    V(i,:) = v_img;
end;
clear v_hdr v_img

% Determine result
%-------------------------------------------------------------------------%
eval(strcat('out = ',expr,';'));