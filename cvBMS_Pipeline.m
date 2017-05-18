% cvBMS Toolkit: Pipeline Template
% _
% Pipeline Template for the cvBMS Toolkit
% For details, see the Manual, page 2-3.
% This script is written for SPM12.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 18/05/2017, 17:45


%%% Step 0: Study parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% project directories
stat_dir = 'C:\Joram\projects\MACS\DataSets\Bogler_et_al_2013\analyses';
work_dir = 'C:\Joram\projects\MACS\DataSets\Bogler_et_al_2013\functional';
% Note: This was the data (Bogler et al., 2013) analyzed in the paper (Soch et al., 2016).

% list of subjects
subj_ids = {'sub01' 'sub02' 'sub03' 'sub04' 'sub05' 'sub06' 'sub07' 'sub08' 'sub09' 'sub10' ...
            'sub11' 'sub12' 'sub13' 'sub14' 'sub15' 'sub16' 'sub17' 'sub18' 'sub19' 'sub20' ...
            'sub21' 'sub22' 'sub23' 'sub24' 'sub25'};

% list of models
mod_nams = {'mod01' 'mod02' 'mod03' 'mod04' 'mod05' 'mod06' 'mod07' 'mod08' 'mod09' 'mod10'};

% model families (optional)
fam_nams = {'fam01' 'fam02' 'fam03' 'fam04'};
mod_fams = [1 1 1 1 2 2 2 3 3 4];

% model space details
ms_name  =  'MS01';
ms_suff  =  'RFX_VB';

% study dimensions
N = numel(subj_ids);
M = numel(mod_nams);
F = numel(fam_nams);


%%% Step 1: First-level model assessment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate model evidences
for i = 1:N                     % all subjects
    for j = 1:M                 % all models
        load(strcat(work_dir,'/',subj_ids{i},'/',mod_nams{j},'/','SPM.mat'));
        MA_cvLME_multi(SPM);
    end;
end;

% get cvLME filename
load(strcat(work_dir,'/',subj_ids{1},'/',mod_nams{1},'/','SPM.mat'));
if numel(SPM.Sess) == 1         % single-session
    data = 'Ky';
    disc = 10 + mod(size(SPM.xX.X,1),10);
    AnC  = false;
    LME_map = strcat('MA_cvLME_',data,'_',int2str(disc),'.nii');
end;
if numel(SPM.Sess) > 1          % multi-session
    data = 'Ky';
    mode = 'N-1';
    AnC  = false;
    LME_map = strcat('MA_cvLME_',data,'_',mode,'.nii');
end;
% Note: This only works if the standard settings from "MA_cvLME_single" and
% "MA_cvLME_multi" are used. If these are modified, the cvLME filename will
% also change.


%%% Step 1.5: Model family inference (optional) %%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate family evidences
for i = 1:N                     % all subjects
    clear LME_maps
    subj_dir = strcat(work_dir,'/',subj_ids{i});
    for j = 1:M                 % all models
        LME_maps{j,1} = strcat(work_dir,'/',subj_ids{i},'/',mod_nams{j},'/',LME_map);
    end;
    MA_LFE_uniform(LME_maps, mod_fams, fam_nams, subj_dir);
end;
% Note: Comment this step, if family inference is not required or desired.


%%% Step 2: Second-level model selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare SPM batch
  BMS_dir = strcat(stat_dir,'/','model_selection','/',ms_name,'_',int2str(M),'mods','_',int2str(N),'subj','_',ms_suff,'/');
% BMS_dir = strcat(stat_dir,'/','model_selection','/',ms_name,'_',int2str(F),'fams','_',int2str(N),'subj','_',ms_suff,'/');
clear sess_map
for i = 1:N                     % all subjects
    clear mod_map
  % clear fam_map
    for j = 1:M % F             % all models
        mod_map{j,1} = strcat(work_dir,'/',subj_ids{i},'/',mod_nams{j},'/',LME_map);
      % fam_map{j,1} = strcat(work_dir,'/',subj_ids{i},'/',fam_nams{j},'/','MA_LFE_uniform.nii');
    end;
    sess_map{1,i}.mod_map = mod_map;
  % sess_map{1,i}.mod_map = fam_map;
end;
% Note: Uncomment the commented lines, if family inference is required.

% create SPM batch
matlabbatch{1}.spm.stats.bms_map.inference.dir{1}      = BMS_dir;
matlabbatch{1}.spm.stats.bms_map.inference.sess_map    = sess_map;
matlabbatch{1}.spm.stats.bms_map.inference.mod_name    = mod_nams';        % fam_nams';
matlabbatch{1}.spm.stats.bms_map.inference.method_maps = 'RFX';            % '<UNDEFINED>';
matlabbatch{1}.spm.stats.bms_map.inference.out_file    = 2;                % 0;
matlabbatch{1}.spm.stats.bms_map.inference.mask        = {''};
matlabbatch{1}.spm.stats.bms_map.inference.nsamp       = '1e6';
filename = strcat(BMS_dir,'design.mat');
save(filename,'matlabbatch');

% perform SPM batch
load(strcat(BMS_dir,'design.mat'));
MS_BMS_group_mods(matlabbatch);
% Note: MS_BMS_group_mods(matlabbatch,[],[],true) invokes
% additional calculation of exceedance probabilities (EPs).


%%% Step 3: Group-level selected-model maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create SMMs
load(strcat(BMS_dir,'BMS.mat'));
MS_SelMod_BMS(BMS);