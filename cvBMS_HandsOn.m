% cvBMS Toolkit: Hands-On Example
% _
% Hands-On Example for the cvBMS Toolkit
% For details, see the Manual, page 4-6.
% This script is written for SPM8.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 15/09/2016, 16:15


%%% Step 0: Pre-processing and model estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Please download the SPM8 Manual or SPM12 Manual and work through Chapter
% 29/31 until the end of Section 29.3/31.3. At this point, you will have
% pre-processed the data and estimated two first-level GLMs, a "categorical
% model" and a "parametric model". [...] After these preliminary analyses,
% the categorical model is located in DIR/categorical and the parametric
% model is located in DIR/parametric where DIR is some folder on your PC.

% analysis directory
DIR = 'C:\Joram\projects\MACS\DataSets\Henson_et_al_2002\stats';


%%% Step 1: First-level model assessment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% categorical model
load(strcat(DIR,'/classical_categorical/SPM.mat'));
MA_cvLME_single(SPM)

% parametric model
load(strcat(DIR,'/classical_parametric/SPM.mat'));
MA_cvLME_single(SPM)


%%% Step 2: "Second-level" model selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bayesian model selection
load(strcat(DIR,'/classical_selection/MS_BMS_group/design.mat'));
MS_BMS_group_mods(matlabbatch)
% Note: MS_BMS_group_mods(matlabbatch,[],[],true) invokes
% additional calculation of exceedance probabilities (EPs).


%%% Step 3: "Group-level" selected-model maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selected-model maps
load(strcat(DIR,'/classical_selection/MS_BMS_group/BMS.mat'));
MS_SelMod_BMS(BMS)