function MS_BMS_group_mods(matlabbatch, method, family, EPs) % Release cvBMS [4]
% _
% Model Inference for Group-Level Bayesian Model Selection
% FORMAT MS_BMS_group_mods(matlabbatch, method, family)
%     matlabbatch - an SPM batch editor job specifying BMS maps
%     method      - a string indicating which methods to use
%     family      - a structure with the following fields:
%     o mods      - a 1 x M vector defining family affiliation (M: models)
%     o fams      - a 1 x F cell array defining family names (F: families)
%     EPs         - a logical indicating exceedance probability calculation
% 
% FORMAT MS_BMS_group_mods(matlabbatch, method, family, EPs) generates
% Bayesian model selection maps according to the batch editor job
% matlabbatch and using the method indicated by method with family
% inference indicated by family.
% 
% The input variable "method" is a string indicating which method to use:
%     If method is 'FFX',    then a fixed effects model is estimated.
%     If method is 'RFX-VB', then a Variational Bayes approach is taken.
%     If method is 'RFX-GS', then a Gibbs Sampling approach is taken.
% 
% The default for this variable is 'RFX-VB' which means that the second-
% level model over log model evidence is estimated using Variational Bayes
% [1,2]. This is recommended if a large space containing lots of voxels is
% analyzed. Setting the method to 'RFX-GS' invokes estimation via Gibbs
% Sampling [3] which is more time-consuming, but also more precise.
% 
% The input variable "family" is an optional structure that specifies
% family inference. The field "mods" is a vector specifying for each model
% which family it belongs to. The field "fams" is a cell array specifying
% the family names. It is required that length(fams) = max(mods). By
% default, this variable is empty.
% 
% When this variable is empty (as by default), a uniform prior over models
% is applied. If this variable is non-empty, this invokes a uniform prior
% over model families in order to allow for unbiased family inference.
% 
% The input variable "EPs" is a logical indicating whether exceedance
% probabilities (EP) are calculated and written to images in case method is
% 'RFX-VB'. If the number of models is very large, it is advisable to only
% calculate EPs at the family level. By default, this variable is false.
% 
% Further information:
%     help ME_BMS_FFX
%     help ME_BMS_RFX_VB
%     help ME_BMS_RFX_GS
%     help MS_BMS_group_fams
%     help MD_Dir_exc_prob
% 
% Exemplary usage:
%     MS_BMS_group_mods(matlabbatch, 'RFX-VB', [], true);
% 
% References:
% [1] Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ (2009):
%     "Bayesian model selection for group studies".
%     NeuroImage, vol. 46, pp. 1004-1017.
% [2] Rosa MJ, Bestmann S, Harrison L, Penny W (2010):
%     "Bayesian model selection maps for group studies".
%     NeuroImage, vol. 49, pp. 217-224.
% [3] Penny WD, Stephan KE, Daunizeau J, Rosa MJ, Friston KJ, Schofield TM,
%     Leff AP (2010): "Comparing Families of Dynamic Causal Models".
%     PLoS ONE, vol. 6, iss. 3, e1000709.
% [4] Soch J, Haynes JD, Allefeld C (2016): "How to avoid mismodelling in
%     GLM-based fMRI data analysis: cross-validated Bayesian model selection".
%     NeuroImage, vol. 141, pp. 469–489.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 05/12/2014, 12:15 (V0.2/V8)
%  Last edit: 15/09/2016, 08:20 (V0.9a/V13a)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get matlabbatch if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    design_mat = spm_select(1,'^*.mat','Select Batch Editor Job!');
    load(design_mat);
    MS_BMS_group_mods(matlabbatch);
    return
else
    if isfield(matlabbatch{1}.spm.stats,'bms')  % catch SPM8 batch version, match SPM12 struct fields
        matlabbatch{1}.spm.stats.bms_map.inference = matlabbatch{1}.spm.stats.bms.bms_map_inf;
        matlabbatch{1}.spm.stats = rmfield(matlabbatch{1}.spm.stats,'bms');
    end;    
    BMS_dir = matlabbatch{1}.spm.stats.bms_map.inference.dir{1};
end;

% Set inference method if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(method), method = 'RFX-VB'; end;

% Set family inference if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(family), family = []; end;

% Inactivate EPs if necessary
%-------------------------------------------------------------------------%
if nargin < 4 || isempty(EPs), EPs = false; end;

% Change to directory if specified
%-------------------------------------------------------------------------%
orig_dir = pwd;
try
    cd(BMS_dir);
catch
    BMS_dir = strcat(pwd,'/');
end

% Get model parameters
%-------------------------------------------------------------------------%
I.sess_map = matlabbatch{1}.spm.stats.bms_map.inference.sess_map;
N = length(I.sess_map);                     % number of subjects
S = length(I.sess_map{1});                  % number of sessions
M = length(I.sess_map{1}(1).mod_map);       % number of models

% Get family dimensions
%-------------------------------------------------------------------------%
if ~isempty(family)
    F  = length(family.fams);
    Mk = sum(repmat([1:F]',[1 M])==repmat(family.mods,[F 1]),2)';
    Mj = Mk(family.mods);
  % F  - the number of families [1 x 1]
  % Mk - the number of models in family k [1 x F]
  % Mj - the number of models in family that model j belongs to [1 x M]
end;

% Get image dimensions
%-------------------------------------------------------------------------%
H = spm_vol(I.sess_map{1}(1).mod_map{1});   % LME image header
V = prod(H.dim);                            % number of voxels

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_BMS_group_mods: load');
spm_progress_bar('Init', 100, 'Load log model evidences...' , '');

% Load log model evidences
%-------------------------------------------------------------------------%
LME = zeros(N,M,V);             % N x M x V array of LMEs
for i = 1:N                     % subjects
    for j = 1:M                 % models
        for k = 1:S             % sessions
            lme_hdr = spm_vol(I.sess_map{i}(k).mod_map{j});
            lme_img = spm_read_vols(lme_hdr);
            lme_img = reshape(lme_img,[1 1 V]);
            % Note: Assuming independence across sessions,
            % log model evidences are summed over sessions.
            LME(i,j,:) = LME(i,j,:) + lme_img;
            spm_progress_bar('Set',(((i-1)*M+(j-1)*S+k)/(N*M*S))*100);
        end;
    end;
end;
clear lme_hdr lme_img

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');

% Create mask image
%-------------------------------------------------------------------------%
LME_1 = squeeze(LME(:,1,:));    % N x V matrix of LMEs for 1st model
if size(LME_1,2) == N, LME_1 = LME_1'; end;
[m_img m_hdr m_ind] = MS_create_mask(LME_1, H);
clear LME_1


%=========================================================================%
% E S T I M A T I O N                                                     %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_BMS_group_mods: estimate');

% Select in-mask voxels only
%-------------------------------------------------------------------------%
v = length(m_ind);
d = floor(v/100);

% Fixed Effects Inference (Bayes' Rule)
%-------------------------------------------------------------------------%
if strcmp(method,'FFX')
    % prior and posterior
    if isempty(family)
        prior = 1/M * ones(1,M);
    else
        prior = 1/F * 1./Mj;
    end;
    PPM   = NaN(M,V);           % posterior probability maps
    % voxel-wise estimation
    spm_progress_bar('Init', 100, 'Estimate posterior probabilities...', '');
    for j = 1:v
        [LGBF, post] = ME_BMS_FFX(LME(:,:,m_ind(j)), prior);
        PPM(:,m_ind(j)) = post';
        if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
    end;
    clear LGBF
end;

% Random Effects Inference (Variational Bayes)
%-------------------------------------------------------------------------%
if strcmp(method,'RFX-VB')
    % prior and posterior
    if isempty(family)
        prior = ones(1,M);
    else
        prior = 1./Mj;
    end;
    alpha = NaN(M,V);           % alpha parameter maps
    LFM   = NaN(M,V);           % likeliest frequency maps
    EFM   = NaN(M,V);           % expected frequency maps
    EPM   = NaN(M,V);           % exceedance probability maps
    % voxel-wise estimation
    spm_progress_bar('Init', 100, 'Estimate posterior distribution...', '');
    for j = 1:v
        alpha(:,m_ind(j)) = ME_BMS_RFX_VB(LME(:,:,m_ind(j)), prior)';
        if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
    end;
    % voxel-wise computation
    spm_progress_bar('Init', 100, 'Compute posterior frequencies...', '');
    for j = 1:v
        LFM(:,m_ind(j)) = MD_Dir_mode(alpha(:,m_ind(j)));
        EFM(:,m_ind(j)) = MD_Dir_mean(alpha(:,m_ind(j)));
        if EPs, EPM(:,m_ind(j)) = MD_Dir_exc_prob(alpha(:,m_ind(j))')'; end;
        if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
    end;
end;

% Random Effects Inference (Gibbs Sampling)
%-------------------------------------------------------------------------%
if strcmp(method,'RFX-GS')
    % prior and posterior
    if isempty(family)
        prior = ones(1,M);
    else
        prior = 1./Mj;
    end;
    alpha = NaN(M,V);           % alpha parameter maps
    EFM   = NaN(M,V);           % expected frequency maps
    EPM   = NaN(M,V);           % exceedance probabilties maps
    % voxel-wise sampling
    spm_progress_bar('Init', 100, 'Sample posterior distribution...', '');
    for j = 1:v
        [alpha_post, exp_freq, exc_prob] = ME_BMS_RFX_GS(LME(:,:,m_ind(j)), prior);
        alpha(:,m_ind(j)) = alpha_post';
        EFM(:,m_ind(j))   = exp_prob';
        EPM(:,m_ind(j))   = exc_prob';
        if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
    end;
end;

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_BMS_group_mods: save');

% Initialise image files
%-------------------------------------------------------------------------%
H = spm_vol(I.sess_map{1}(1).mod_map{1});

% Retrieve model names
%-------------------------------------------------------------------------%
mods = matlabbatch{1}.spm.stats.bms_map.inference.mod_name';

% Write images to disk
%-------------------------------------------------------------------------%
for i = 1:M
    if strcmp(method,'FFX')     % Bayes' Rule
        H.fname   = strcat(mods{i},'_model_PPM.nii');
        H.descrip = 'MS_BMS_group_mods: posterior probability maps';
        spm_write_vol(H,reshape(PPM(i,:),H.dim));
        BMS.map.ffx.ppm{i} = H.fname;
    end;
    if strcmp(method,'RFX-VB')  % Variational Bayes
        H.fname   = strcat(mods{i},'_model_alpha.nii');
        H.descrip = 'MS_BMS_group_mods: alpha parameters';
        spm_write_vol(H,reshape(alpha(i,:),H.dim));
        BMS.map.rfx.alpha{i} = H.fname;
        H.fname   = strcat(mods{i},'_model_LFM.nii');
        H.descrip = 'MS_BMS_group_mods: likeliest frequencies';
        spm_write_vol(H,reshape(LFM(i,:),H.dim));
        BMS.map.rfx.lfm{i} = H.fname;
        H.fname   = strcat(mods{i},'_model_EFM.nii');
        H.descrip = 'MS_BMS_group_mods: expected frequencies';
        spm_write_vol(H,reshape(EFM(i,:),H.dim));
        BMS.map.rfx.efm{i} = H.fname;
        if EPs
            H.fname   = strcat(mods{i},'_model_EPM.nii');
            H.descrip = 'MS_BMS_group_mods: exceedance probabilities';
            spm_write_vol(H,reshape(EPM(i,:),H.dim));
            BMS.map.rfx.epm{i} = H.fname;
        end;
    end;
    if strcmp(method,'RFX-GS')  % Gibbs Sampling
        H.fname   = strcat(mods{i},'_model_alpha.nii');
        H.descrip = 'MS_BMS_group_mods: alpha parameters';
        spm_write_vol(H,reshape(alpha(i,:),H.dim));
        BMS.map.rfx.alpha{i} = H.fname;
        H.fname   = strcat(mods{i},'_model_EFM.nii');
        H.descrip = 'MS_BMS_group_mods: expected frequencies';
        spm_write_vol(H,reshape(EFM(i,:),H.dim));
        BMS.map.rfx.efm{i} = H.fname;
        H.fname   = strcat(mods{i},'_model_EPM.nii');
        H.descrip = 'MS_BMS_group_mods: exceedance probabilities';
        spm_write_vol(H,reshape(EPM(i,:),H.dim));
        BMS.map.rfx.epm{i} = H.fname;
    end;
end;

% Write mask image to disk
%-------------------------------------------------------------------------%
spm_write_vol(m_hdr,reshape(m_img,H.dim));

% Finalize BMS structure
%-------------------------------------------------------------------------%
if strcmp(method,'FFX')
    BMS.map.ffx.data  = matlabbatch{1}.spm.stats.bms_map.inference.sess_map;
    BMS.map.ffx.mask  = m_hdr.fname;
    BMS.map.ffx.prior = prior;
end;
if strcmp(method,'RFX-VB') || strcmp(method,'RFX-GS')
    BMS.map.rfx.data  = matlabbatch{1}.spm.stats.bms_map.inference.sess_map;
    BMS.map.rfx.mask  = m_hdr.fname;
    BMS.map.rfx.prior = prior;
end;
BMS.fname = strcat(BMS_dir,'BMS.mat');
save(strcat(BMS_dir,'BMS.mat'),'BMS');

% Perform family inference
%-------------------------------------------------------------------------%
if ~isempty(family)
    MS_BMS_group_fams(BMS, family.mods, family.fams);
end;

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);