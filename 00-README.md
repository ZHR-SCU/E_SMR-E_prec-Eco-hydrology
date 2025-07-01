# E_SMR-E_prec-Eco-hydro
################################################################################## 
##################################################################################
# This research code is mainly divided into MATLAB code and Python code. We used MATLAB for data processing, analysis of θcri, implementation of sliding segmentation, and modeling of non-linear mixed effects model in the main research. We use Python to quantify the contribution rate of external moisture meteorological element variables to soil moisture content at different levels.
##################################################################################
      
      Matlab code section

            This section outlines the overall workflow for the MATLAB R2022b scripts used in this study.

                  Step 1: Identify the site-specific θ_cri (critical soil moisture threshold)
                     using the procedures provided in "1 - theta_cri_determination.m".
 
                Step 2: Visualize the soil moisture data and perform manual outlier checks using
                          "2 - extractAndVisualizeSMData.m". This script also highlights
                                       periods when soil moisture exceeds θ_cri.
 
                Step 3: Based on the outcomes from Step 2, apply "3 - analyzeESmrprecEffect.m" 
                    to perform a segmented analysis of spring soil moisture enhancement,
                                    and quantify both E_SMR and E_prec.
 
                Step 4: Evaluate the impacts of E_SMR, E_prec, and near-surface air temperature (airT) 
                    on the start of season (SOS) using the nonlinear mixed-effects model, 
                                and quantify their relative contributions.
          
                       We uploaded a Demon.mat database to test the usability of the code
           
                                   All the code was written by Hanrui Zhao. 
                    For any questions, please contact the author at hanrui_zhao@stu.scu.edu.cn
################################################################################## 

      Python code section
                      This project demonstrates how to use SHAP (SHapley Additive exPlanations) to interpret
                        the influence of hydrometeorological variables, such as precipitation and active layer thickness (ALT), 
                        on different layers of soil moisture using a Random Forest regression model.
                        
                      The python code demonstrates how to use SHAP (SHapley Additive exPlanations) to interpret the influence of hydrometeorological variables, 
                         such as precipitation and active layer thickness (ALT), on different layers of soil moisture using a Random Forest regression model.

                      We trained two Random Forest regression models to predict:
                                    - Shallow soil moisture (SM_bottom)
                                    - Root-zone soil moisture (SM_root)
﻿
                      Each model takes the following predictors:
                                    - Standardized precipitation (Prec)
                                    - Standardized active layer thickness (ALT)
                                    - Interaction term of precipitation and ALT
﻿
                       These features were preprocessed using z-score normalization, and the response variables were also standardized.
﻿
                       Rather than splitting the data into training and test sets, we trained the models on the full available dataset.
                       This approach was adopted because our objective was not to optimize prediction accuracy, but to perform explanatory modeling, specifically, to                                 quantify the relative contribution of each predictor to soil moisture variability using SHAP analysis.
﻿
                       Training the model on the full dataset provides a more stable and comprehensive interpretation of variable influence,
                       which is particularly suitable for limited or observational datasets where test-train splits may reduce the robustness of inference.
﻿
################################################################################## 
