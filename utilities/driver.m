 PATH_TO_PYTHON = '/Users/Kaiming Xu/anaconda3/envs/rOMT';
addpath('./Sensitivities', './analyzeFlows', genpath('./utilities'))
% Select dataset (ONE OR MULTIPLE)

TAG = 'C294_small'; % brain (Chen et al. 2022)

% %% RUN ONE DATASET
 [experiment, flag] = runRegularizedOMT(TAG);
% process and visualize results
runPostprocessing(TAG, PATH_TO_PYTHON);




