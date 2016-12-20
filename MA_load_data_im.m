function Y = MA_load_data_im(SPM, m_ind)
% _
% Load Data from General Linear Model (in-mask)
% FORMAT Y = MA_load_mask_im(SPM, m_ind)
% 
%     SPM   - a structure specifying an estimated GLM
%     m_ind - a  1 x v vector indexing in-mask voxels
% 
%     Y     - an n x v data matrix (n: scans; v: voxels)
% 
% FORMAT Y = MA_load_data_im(SPM, m_ind) loads only in-mask time series
% belonging to the GLM specified by SPM and returns a data matrix.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 24/10/2014, 18:30 (V0.2/V6)
%  Last edit: 19/11/2015, 05:00 (V0.4/V12)


% Get data dimensions
%-------------------------------------------------------------------------%
n = length(SPM.xY.VY);
v = numel(m_ind);
d = ceil(n/100);

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_load_data_im: load');
spm_progress_bar('Init',100,'Load in-mask time series...','');

% Load time series
%-------------------------------------------------------------------------%
Y = zeros(n,v);
for i = 1:n
    y_img  = spm_read_vols(SPM.xY.VY(i));
  % y_img  = reshape(y_img,[1 prod(SPM.VM.dim)]);
    Y(i,:) = y_img(m_ind);
    if mod(i,d) == 0, spm_progress_bar('Set',(i/n)*100); end;
end;

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');