%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This matlab code is used to construct the nonlinear mixed-effects (NLME) modelconsidering regional differences in this study, 
%     including aspects such as the standardization of input predictor variables, the modeling, model diagnosis, 
%                 data orthologization, contribution rate analysis, and plotting.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all

%% Data preparation
E_SMR = E_SMR;
E_prec = E_prec;
SOS = SOS;
airT = airT;

%% Modeling
% Suppose there are three regions, observed for 10 years, 12 years and 9 years respectively
n1 = 10;
n2 = 12;
n3 = 9;
regionID = [ repelem("S1",n1).'; ...
           repelem("S2",n2).'; ...
           repelem("S3",n3).'; ...
           ];

% Make it categorical
regionCat = categorical(regionID);

% Data standardization
E_SMR = zscore(E_SMR, 0, 'omitnan');
E_prec = zscore(E_prec, 0, 'omitnan');
airT = zscore(airT, 0, 'omitnan');

% Bulding table 
tbl = table( ...
    SOS, ...
    E_SMR, ...
    E_prec, ...
    airT, ...
    'VariableNames', {'SOS','E_SMR','E_prec','airT'} ...
    );
tbl.region = regionCat;

% modeling
formula = 'SOS ~ 1 + (E_SMR + E_prec + airT)^3+(1|region)';
lme = fitlme(tbl, formula);
disp(lme); % Check results

% Visualization of NLME marginal results (Taking predictor E_SMR as an example)
airT_levels  = [-1, 0, 1];
E_prec_levels = [-1, 0, 1];
E_SMR_levels = [-1, 0, 1];
E_SMR_vals = linspace(min(dataTable_clean.E_SMR), max(dataTable_clean.E_SMR), 50);
figure;
colors = lines(length(airT_levels)); 
line_styles = {'-', '--', ':'}; 
legend_entries = {};
hold on;
idx = 1;
for i = 1:length(airT_levels)
    for j = 1:length(E_prec_levels)
        airT_val = airT_levels(i);
        E_prec_val = E_prec_levels(j);      
        newData_int = table();
        newData_int.E_SMR = E_SMR_vals(:);
        newData_int.E_prec = repmat(E_prec_val, length(E_SMR_vals), 1);
        newData_int.airT = repmat(airT_val, length(E_SMR_vals), 1);
        newData_int.Elevation = repmat(Elevation_fixed, length(E_SMR_vals), 1);
        % Since the marginal effect does not consider regional differences, 
        % a virtual site is added (because there is (1|region) in the model, this column must be present).
        newData_int.site = repmat(categorical("dummy_site"), length(E_SMR_vals), 1);       
        % Predicted SOS
        SOS_int = predict(lme, newData_int);        
        % Draw the curve (the color represents airT and the line type represents E_prec)
        plot(E_SMR_vals, SOS_int, 'Color', colors(i,:), 'LineWidth', 2, 'LineStyle', line_styles{j});        
        % Update the legend
        legend_entries{idx} = sprintf('{\it{T}}_{air} = %.1f, {\it{E}}_{SMR} = %.1f', airT_val, E_prec_val);
        idx = idx + 1;
    end
end
xlabel('{\it{E}}_{SMR} (-)');
ylabel('Predicted SOS (DOY)');
legend(legend_entries, 'Location', 'bestoutside');
grid on;
hold off;

%% Model diagnosis
residualss = lme.Residuals.Raw;
% Q-Q plot 
figure;
qqplot(residualss);
% Residual and fitted value scatter plot
figure;
plotResiduals(lme, 'fitted');
% Cook distance
res = residuals(lme, 'ResidualType', 'Standardized'); % Extract necessary variables
X = lme.designMatrix('Fixed'); 
n = size(X, 1);                         % Number of observations
p = size(X, 2);                         % Number of fixed effect parameters
% Diagonal elements of the hat matrix (leverage values)
H = X * ((X' * X) \ X');               % Hat matrix
h = diag(H);                           % Leverage values
% Residual sum of squares
MSE = sum(res.^2) / (n - p);           % Mean squared error
% Compute Cook's distance
cookD = (res.^2 ./ (p * MSE)) .* (h ./ ((1 - h).^2));
% Plot Cook's distance
figure;
stem(cookD, 'filled');
% Mark common threshold
hold on;
threshold = 4 / n;
yline(threshold, 'r--');
title("Cook's Distance");
xlabel('Observation');
ylabel("Cook's Distance");
legend('Cook''s D', 'Threshold = 4/n');

%% Check the marginal and conditional R^2
[psi, mse] = covarianceParameters(lme);
Sigma_rand = psi{1};           % Covariance matrix of random intercept (1×1)
sigma_rand = Sigma_rand(1,1);  % Variance of random effects σ^2_rand
sigma_res  = mse;              % Residual variance σ^2_res

% Compute the variance of fitted values from fixed effects
% Extract fitted values based only on fixed effects (exclude random effects)
yhat_fixed = fitted(lme, 'Conditional', false);
%  Compute the total variance (using N in the denominator)
sigma_fix  = var(yhat_fixed, 1);  % σ^2_fix

% Calculate R^2 following Nakagawa & Schielzeth (2013)
R2_marginal    = sigma_fix / (sigma_fix + sigma_rand + sigma_res);
R2_conditional = (sigma_fix + sigma_rand) / (sigma_fix + sigma_rand + sigma_res);
fprintf('Marginal R^2 (fixed effects only) = %.3f\n', R2_marginal);
fprintf('Conditional R^2 (fixed + random) = %.3f\n', R2_conditional);

%% Visualization of conditional prediction results of NLME model (Each region)
LineStyles = {'-', '--', ':'};
colors = lines(length(airT_levels));
region_categories = categorical({'S1','S2','S3'}, true);
figure
for s = 1:length(region_categories)
    this_region = region_categories(s); 
    subplot(rows, cols, s);
    hold on;
    grid on;
    box on;
    for i = 1:length(airT_levels)
        for j = 1:length(E_SMR_levels)
            newData = table;
            newData.E_SMR = E_SMR_vals(:);
            newData.E_prec = repmat(E_SMR_levels(j), length(E_prec_vals), 1);
            newData.airT = repmat(airT_levels(i), length(E_SMR_vals), 1);
            newData.Elevation = repmat(Elevation_fixed, length(E_SMR_vals), 1);
            newData.site = repmat(categorical(this_region, region_categories), height(newData), 1);
            SOS_pred = predict(lme, newData);
            plot(E_SMR_vals, SOS_pred, ...
                'Color', colors(i,:), ...
                'LineWidth', 1.5, ...
                'LineStyle', LineStyles{j} ...
                );
        end
    end
    title(char(this_region), 'FontSize', 9);
end
han = axes(gcf, 'visible', 'off');
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
xlabel(han, '{\it{E}}_{SMR} (-)');
ylabel(han, 'Predicted SOS (DOY)');

%% Data orthogonalization and contribution quantification
X_raw = [ ...
    E_SMR,              E_prec,             airT,...                  % degree 1
    E_SMR.^2,           E_SMR.*E_prec,       E_SMR.*airT, ...           % degree 2
    E_prec.^2,          E_prec.*airT,       airT.^2, ...               % degree 2
    E_SMR.^3,           E_SMR.^2.*E_prec,    E_SMR.^2.*airT, ...        % degree 3
    E_SMR.*E_prec.^2,    E_SMR.*E_prec.*airT, E_SMR.*airT.^2, ...        % degree 3
    E_prec.^3,          E_prec.^2.*airT,    E_prec.*airT.^2, ...        % degree 3
    airT.^3                                           ...             % degree 3
];
% Perform QR decomposition on X_raw to obtain orthogonal basis Q
[Q, ~] = qr(X_raw, 0);    % Q(:,k)' * Q(:,j) = 0 (j ≠ k), and norm(Q(:,k)) = 1

% Construct a table for model fitting
% Assign names Q1, Q2, ..., Qn to columns of Q
nQ = size(Q, 2);
Qnames = arrayfun(@(k) sprintf('Q%d', k), 1:nQ, 'UniformOutput', false);

% Combine SOS and Q into one table
tbl = array2table([SOS, Q], 'VariableNames', ['SOS', Qnames]);

% Fit the linear model using orthogonalized predictors
% Q1...Qn already include all higher-order and interaction terms
formulaQ = ['SOS ~ 1 + ', strjoin(Qnames, ' + ')];
mdl = fitlm(tbl, formulaQ);

% Check variance inflation factors (VIF) to confirm orthogonalization effectiveness
disp(mdl);
anovaTbl = anova(mdl, 'summary');
disp(anovaTbl);
X = mdl.Variables(:, 2:end);         % Exclude response variable
X = table2array(X);                  % Convert to numeric matrix
VIF = zeros(1, size(X, 2));          % Initialize VIF array
for i = 1:size(X, 2)
    Xi = X(:, i);
    X_other = X(:, [1:i-1, i+1:end]);
    mdl_i = fitlm(X_other, Xi);
    R2_i = mdl_i.Rsquared.Ordinary;
    VIF(i) = 1 / (1 - R2_i);         % Calculate VIF
end
disp('Variance Inflation Factors (VIF):');
disp(array2table(VIF, 'VariableNames', Qnames));

% Calculate relative contribution of each orthogonalized predictor
coeff = mdl.Coefficients.Estimate(2:end);  % Exclude intercept
SS_total = sum((mdl.Fitted - mean(mdl.Fitted)).^2);
contribution = (coeff.^2) ./ sum(coeff.^2) * mdl.Rsquared.Ordinary;

% Display contribution table
disp('Relative Contribution of Each Orthogonal Predictor:');
contribution_table = table(Qnames', contribution', 'VariableNames', {'Predictor', 'Contribution'});
disp(contribution_table);





