function [T_AUC, T_ACC, T_PPV, T_NPV] = feature_MOR_AUC_lassoselection(X_baseline,X_romics,fnames_romics,fnames_baseline,y_target)
% Here, build the feature versus log MOR and ROC AUC tables. 
% Also has feature selection frequency from bootstrapping.

% Model A is the baseline. 
% Model B is the radiomic model. 
% Model C is the combined baseline+radiomic model
%% load the data to work with
% % % % % % disp('Loading data...')
% % % % % % load(['C:\Users\zakar\OneDrive\Documents\MATLAB\Radiomics\Data_extraction_bin\X_y_data_' num2str(num) '.mat']);
% % % % % % 
% % % % % % eval(['X_baseline = X_baseline_T2_' num2str(num) ';']);
% % % % % % eval(['X_romics   = X_romics_T2_' num2str(num) ';']);
% % % % % % eval(['fnames_romics  = fnames_romics_' num2str(num) ';']);
% % % % % % eval(['fnames_baseline  = fnames_baseline_' num2str(num) ';']);
% % % % % % eval(['y_target   = y_target_' num2str(num) ';']);

%% model A: baseline only
disp('Baseline analysis starting...')

X = [X_baseline];

% what are the baseline feature names?
fnames = [fnames_baseline];
y = y_target;

mdl_A = fitclinear(X,y,'Learner','logistic','Regularization','lasso');

[label,scores] = predict(mdl_A,X);
[Xroc,Yroc,~,AUC_A] = perfcurve(y, scores(:,2), 1);

beta_coefficients_fnames_A = [fnames',num2cell(mdl_A.Beta)];

% k-fold cross validation, k=5
[AUC_all_lasso_A, Accuracy_all_A, PPV_all_A, NPV_all_A, fnames_beta_coefficients_kfold_A, mdl_median_fold_A, scores_median_fold_A] = ...
    linear_model_validation(X, y, fnames, 5, 'lasso');

% bootstrap part
% AUC_lasso_bootstrap_A = bootstrp(100, @(X) {get_roc_auc(X, y, 'lasso')}, X);

% AUC_avg_lasso_bootstrapdata_A = ...
%     bootstrp(100, @(X) {linear_model_validation_average(X, y, 5, 'lasso')}, X);

%% model B: radiomics model only
disp('Radiomics analysis starting...')

X = X_romics;
fnames = fnames_romics;
y = y_target;

mdl_B = fitclinear(X,y,'Learner','logistic','Regularization','lasso');

[label,scores] = predict(mdl_B,X);
[Xroc,Yroc,~,AUC_B] = perfcurve(y, scores(:,2), 1);

beta_coefficients_fnames_B = [fnames',num2cell(mdl_B.Beta)];

% k-fold cross validation, k=5
[AUC_all_lasso_B, Accuracy_all_B, PPV_all_B, NPV_all_B, fnames_beta_coefficients_kfold_B, mdl_median_fold_B, scores_median_fold_B] = linear_model_validation(X, y, fnames, 5, 'lasso');

% bootstrap part
% AUC_lasso_bootstrap_B = bootstrp(100, @(X) {get_roc_auc(X, y, 'lasso')}, X);

% AUC_avg_lasso_bootstrapdata_B = ...
%     bootstrp(100, @(X) {linear_model_validation_average(X, y, 5, 'lasso')}, X);

%% model C: combined model. 
disp('Combined model starting...')

X = [X_romics, X_baseline];
fnames = [fnames_romics, fnames_baseline];
y = y_target;

mdl_C = fitclinear(X,y,'Learner','logistic','Regularization','lasso');

[label,scores] = predict(mdl_C,X);
[Xroc,Yroc,~,AUC_C] = perfcurve(y, scores(:,2), 1);

beta_coefficients_fnames_C = [fnames',num2cell(mdl_C.Beta)];

% k-fold cross validation, k=5
[AUC_all_lasso_C, Accuracy_all_C, PPV_all_C, NPV_all_C, fnames_beta_coefficients_kfold_C, mdl_median_fold_C, scores_median_fold_C] = linear_model_validation(X, y, fnames, 5, 'lasso');

% bootstrap part
%AUC_lasso_bootstrap_C = bootstrp(100, @(X) {get_roc_auc(X, y, 'lasso')}, X);

%AUC_avg_lasso_bootstrapdata_C = ...
    %bootstrp(100, @(X) {linear_model_validation_average(X, y, 5, 'lasso')}, X);

% Save the model and scores. 
save('mdl_scores_LASSO.mat','mdl_median_fold_A','scores_median_fold_A','mdl_median_fold_B','scores_median_fold_B','mdl_median_fold_C','scores_median_fold_C');

%% More validation
% % % % % % % % % % % % Using bootstrapping to see feature selection frequency in the combined
% % % % % % % % % % % % model.
% % % % % % % % % % % bootstrapdata_lasso = bootstrp(1000, @(X) {lasso_selection_frequency(X, y, fnames, 5)}, X);
% % % % % % % % % % % whole_array = [];
% % % % % % % % % % % for i = 1:length(bootstrapdata_lasso),
% % % % % % % % % % %     whole_array = [whole_array; bootstrapdata_lasso{i,1}];
% % % % % % % % % % % end
% % % % % % % % % % % x_hist = fnames;
% % % % % % % % % % % y_hist = sum(whole_array);
% % % % % % % % % % % [y_hist_descend, sort_idx] = sort(y_hist,'descend');
% % % % % % % % % % % fnames = fnames(sort_idx);
% % % % % % % % % % % 
% % % % % % % % % % % figure;
% % % % % % % % % % % X = categorical(fnames);
% % % % % % % % % % % X = reordercats(X,fnames);
% % % % % % % % % % % bar(X,y_hist_descend,0.2)
% % % % % % % % % % % %title([' Cross-validated expected frequency of LASSO feature selection']);
% % % % % % % % % % % ylabel('Frequency');
% % % % % % % % % % % set(get(gca,'YLabel'),'Rotation',90)
% % % % % % % % % % % ytickangle(90);

%% save all the tables
% feature versus log-MOR.
% you'll have fname -  all data logMOR(which is beta) - kfold logMOR (beta).
flogMORtable_T2_LASSO = [[{'/'},{'All'},{'1'},{'2'},{'3'},{'4'},{'5'}];
                   [{'Model A'},{'/'},{'/'},{'/'},{'/'},{'/'},{'/'}];
                   [beta_coefficients_fnames_A,fnames_beta_coefficients_kfold_A(:,2:end)];
                   [{'Model B'},{'/'},{'/'},{'/'},{'/'},{'/'},{'/'}];
                   [beta_coefficients_fnames_B,fnames_beta_coefficients_kfold_B(:,2:end)];
                   [{'Model C'},{'/'},{'/'},{'/'},{'/'},{'/'},{'/'}];
                   [beta_coefficients_fnames_C,fnames_beta_coefficients_kfold_C(:,2:end)]];
save('flogMORdata_LASSO.mat','flogMORtable_T2_LASSO');

%% summary table
disp('Writing table...')
T_AUC={mean(AUC_A) [] mean(AUC_B) [] mean(AUC_C) []};
% trained on all data expected AUC
% T_AUC=[T_AUC;{mean(cell2mat(AUC_lasso_bootstrap_A)) std(cell2mat(AUC_lasso_bootstrap_A)) mean(cell2mat(AUC_lasso_bootstrap_B)) std(cell2mat(AUC_lasso_bootstrap_B)) mean(cell2mat(AUC_lasso_bootstrap_C)) std(cell2mat(AUC_lasso_bootstrap_C))}];
% cross-fold training average AUC
T_AUC=[T_AUC;{mean(AUC_all_lasso_A) std(AUC_all_lasso_A) mean(AUC_all_lasso_B) std(AUC_all_lasso_B) mean(AUC_all_lasso_C) std(AUC_all_lasso_C)}];
% cross-fold training expected average AUC
% T_AUC=[T_AUC;{mean(cell2mat(AUC_avg_lasso_bootstrapdata_A)) std(cell2mat(AUC_avg_lasso_bootstrapdata_A)) mean(cell2mat(AUC_avg_lasso_bootstrapdata_B)) std(cell2mat(AUC_avg_lasso_bootstrapdata_B)) mean(cell2mat(AUC_avg_lasso_bootstrapdata_C)) std(cell2mat(AUC_avg_lasso_bootstrapdata_C))}];
disp('Done.')
%
T_ACC=[];
T_ACC=[{'Mean_Accuracy_Model_A:'},{'Stdev_Accuracy_Model_A:'}];
T_ACC=[T_ACC;{mean(Accuracy_all_A) std(Accuracy_all_A)}];
T_ACC=[T_ACC;[{'Mean_Accuracy_Model_B'},{'Stdev_Accuracy_Model_B'}]];
T_ACC=[T_ACC;{mean(Accuracy_all_B) std(Accuracy_all_B)}];
T_ACC=[T_ACC;[{'Mean_Accuracy_Model_C'},{'Stdev_Accuracy_Model_C'}]];
T_ACC=[T_ACC;{mean(Accuracy_all_C) std(Accuracy_all_C)}];
%
T_PPV=[];
T_PPV=[{'Mean_PPV_Model_A:'},{'Stdev_PPV_Model_A:'}];
T_PPV=[T_PPV;{mean(PPV_all_A) std(PPV_all_A)}];
T_PPV=[T_PPV;[{'Mean_PPV_Model_B'},{'Stdev_PPV_Model_B'}]];
T_PPV=[T_PPV;{mean(PPV_all_B) std(PPV_all_B)}];
T_PPV=[T_PPV;[{'Mean_PPV_Model_C'},{'Stdev_PPV_Model_C'}]];
T_PPV=[T_PPV;{mean(PPV_all_C) std(PPV_all_C)}];
%
T_NPV=[];
T_NPV=[{'Mean_NPV_Model_A:'},{'Stdev_NPV_Model_A:'}];
T_NPV=[T_NPV;{mean(NPV_all_A) std(NPV_all_A)}];
T_NPV=[T_NPV;[{'Mean_NPV_Model_B'},{'Stdev_NPV_Model_B'}]];
T_NPV=[T_NPV;{mean(NPV_all_B) std(NPV_all_B)}];
T_NPV=[T_NPV;[{'Mean_NPV_Model_C'},{'Stdev_NPV_Model_C'}]];
T_NPV=[T_NPV;{mean(NPV_all_C) std(NPV_all_C)}];
end