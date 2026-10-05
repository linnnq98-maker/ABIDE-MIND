%% =========================================================
% BrainSMASH-based spatial permutation test
% Association between PLS1 regional scores and imaging map
% =========================================================

% --------------------------
% Parameter settings
% --------------------------
nregs_lh = 152;   % Number of left-hemisphere regions


% --------------------------
% Load data
% --------------------------

% PLS1 regional scores
PLS1_scores = readmatrix( ...
    fullfile('path', 'to', 'data', 'PLS1_ROIscores.csv') ...
);

% Retain left-hemisphere regions
PLS1_scores_lh = PLS1_scores(1:nregs_lh);


% Observed imaging-statistic map
Y_true = readmatrix( ...
    fullfile('path', 'to', 'data', 'subtype1_vs_control_glm_t_value.csv') ...
);

% Retain left-hemisphere regions
Y_true_lh = Y_true(1:nregs_lh);


% BrainSMASH surrogate maps
% Expected format: N permutations x number of regions
surrogates = readmatrix( ...
    fullfile('path', 'to', 'data', 'subtype1_surrogate.csv') ...
);

% Retain left-hemisphere regions
surrogates_lh = surrogates(:, 1:nregs_lh);


% --------------------------
% Check input dimensions
% --------------------------
assert(length(PLS1_scores) >= nregs_lh, ...
    'PLS1 score file contains fewer than 152 regions.');

assert(length(Y_true) >= nregs_lh, ...
    'Observed imaging map contains fewer than 152 regions.');

assert(size(surrogates, 2) >= nregs_lh, ...
    'Surrogate maps contain fewer than 152 regions.');


% --------------------------
% Z-score normalization
% --------------------------
PLS1_scores_lh = zscore(PLS1_scores_lh);

Y_true_lh = zscore(Y_true_lh);

% Normalize each surrogate map separately
surrogates_lh = zscore(surrogates_lh, 0, 2);


% --------------------------
% Observed correlation
% --------------------------
R_real = corr(PLS1_scores_lh, Y_true_lh);


% --------------------------
% Surrogate correlation distribution
% --------------------------
n_perm = size(surrogates_lh, 1);

R_perm = zeros(n_perm, 1);

for i = 1:n_perm
    
    R_perm(i) = corr( ...
        PLS1_scores_lh, ...
        surrogates_lh(i, :)' ...
    );
    
end


% --------------------------
% Two-sided permutation P value
% --------------------------
% +1 correction prevents an empirical P value of zero
p_perm = ...
    (sum(abs(R_perm) >= abs(R_real)) + 1) / ...
    (n_perm + 1);


% --------------------------
% Save results
% --------------------------
writematrix( ...
    R_real, ...
    'subtype1_real_corr.csv' ...
);

writematrix( ...
    R_perm, ...
    'subtype1_permutation_corrs.csv' ...
);

result_table = table( ...
    R_real, ...
    p_perm, ...
    n_perm, ...
    'VariableNames', ...
    {'Observed_r', 'Permutation_p', 'N_permutations'} ...
);

writetable( ...
    result_table, ...
    'subtype1_brainsmash_permutation_result.csv' ...
);


% --------------------------
% Print results
% --------------------------
fprintf( ...
    ['PLS1 correlation with the observed imaging map: ', ...
     'r = %.4f, permutation p = %.5f\n'], ...
    R_real, ...
    p_perm ...
);


% --------------------------
% Visualize null distribution
% --------------------------
figure;

histogram( ...
    R_perm, ...
    'Normalization', 'pdf', ...
    'FaceColor', [0.5 0.5 0.5] ...
);

hold on;

line( ...
    [R_real R_real], ...
    ylim, ...
    'Color', 'r', ...
    'LineWidth', 2 ...
);

xlabel('Correlation coefficient');

ylabel('Probability density');

title( ...
    sprintf( ...
        'BrainSMASH permutation test: r = %.3f, p = %.5f', ...
        R_real, ...
        p_perm ...
    ) ...
);

set(gca, 'FontSize', 12);

grid on;