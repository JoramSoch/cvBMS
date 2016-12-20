function Y = MA_load_data(SPM)
% _
% Load Data from General Linear Model
% FORMAT Y = MA_load_mask(SPM)
% 
%     SPM   - a structure specifying an estimated GLM
% 
%     Y     - an n x v data matrix (n: scans; v: voxels)
% 
% FORMAT Y = MA_load_data(SPM) loads the time series belonging to the GLM
% specified by SPM and returns a data matrix.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 24/10/2014, 18:30 (V0.2/V6)
%  Last edit: 27/11/2014, 17:00 (V0.2/V8)


% Get data dimensions
%-------------------------------------------------------------------------%
n = length(SPM.xY.VY);
v = prod(SPM.VM.dim);
d = ceil(n/100);

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_load_data: load');
spm_progress_bar('Init',100,'Load time series...','');

% Load time series
%-------------------------------------------------------------------------%
Y = zeros(n,v);
for i = 1:n
    y_img  = spm_read_vols(SPM.xY.VY(i));
    y_img  = reshape(y_img,[1 v]);
    Y(i,:) = y_img;
    if mod(i,d) == 0, spm_progress_bar('Set',(i/n)*100); end;
end;

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');