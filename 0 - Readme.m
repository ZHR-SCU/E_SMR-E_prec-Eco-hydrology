%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% 
% This section outlines the overall workflow for the MATLAB scripts used in this study.
% 
%    Step 1: Identify the site-specific θ_cri (critical soil moisture threshold)
%         using the procedures provided in "1 - theta_cri_determination.m".
% 
%    Step 2: Visualize the soil moisture data and perform manual outlier checks using
%        "2 - extractAndVisualizeSMData.m". This script also highlights
%                     periods when soil moisture exceeds θ_cri.
% 
%    Step 3: Based on the outcomes from Step 2, apply "3 - analyzeESmrprecEffect.m" 
%        to perform a segmented analysis of spring soil moisture enhancement,
%                        and quantify both E_SMR and E_prec.
% 
%    Step 4: Evaluate the impacts of E_SMR, E_prec, and near-surface air temperature (airT) 
%        on the start of season (SOS) using the nonlinear mixed-effects model, 
%                    and quantify their relative contributions.
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       All the code was written by Hanrui Zhao. 
%        For any questions, please contact the author at hanrui_zhao@stu.scu.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
