%% Unified MATLAB Simulation Framework for Atomic Linearized Projections
clc; clear; close all;

%% 1. Load Dataset from Excel
file_name = 'Atomic Data.xlsx'; 
if exist(file_name, 'file') ~= 2
    error('Error: "%s" not found! Please ensure it is in the active folder.', file_name);
end

data = readtable(file_name);
Z = data{:, 1}; % Column 1: Atomic Number (Z)
A = data{:, 2}; % Column 2: Atomic Mass (A)

%% 2. Core Formulations & Parametric Boundary (L2)
h = ((A - Z) ./ (A + Z)).^2;
L2 = 2 * sqrt(4 * (A.^2 - Z.^2) + pi^2 * Z.^2);

%% 3. Compute 8 Heuristic Structural Combinations (Cases I - VIII)
base_R = pi * (A + Z);
Epsilon_Matrix = zeros(length(Z), 8);

% Cases 1 to 4 (To be plotted against Atomic Number - Z)
Epsilon_Matrix(:, 1) = base_R .* (1 + (3*h)./(10 + sqrt(4 - 3*h))) - L2; 
Epsilon_Matrix(:, 2) = base_R .* (1 + (3*h)./(10 - sqrt(4 + 3*h))) - L2; 
Epsilon_Matrix(:, 3) = base_R .* (1 + (3*h)./(10 + sqrt(4 + 3*h))) - L2; 
Epsilon_Matrix(:, 4) = base_R .* (1 - (3*h)./(10 + sqrt(4 - 3*h))) - L2; 

% Cases 5 to 8 (To be plotted against Atomic Mass - A)
Epsilon_Matrix(:, 5) = base_R .* (1 - (3*h)./(10 - sqrt(4 + 3*h))) - L2; 
Epsilon_Matrix(:, 6) = base_R .* (1 + (3*h)./(10 + sqrt(4 + 3*h))) - L2; 
Epsilon_Matrix(:, 7) = base_R .* (1 + (3*h)./(10 - sqrt(4 - 3*h))) - L2; % [OPTIMIZED LINEAR CASE]
Epsilon_Matrix(:, 8) = base_R .* (1 - (3*h)./(10 - sqrt(4 - 3*h))) - L2; 

%% 4. Statistical Evaluation & Linear Fitting Metrics
fprintf('\n========================================================================\n');
fprintf('                STATISTICAL METRICS FOR ALL HEURISTIC CASES             \n');
fprintf('========================================================================\n');
fprintf('%-10s %-12s %-10s %-10s %-10s\n', 'Case', 'Pearson_R', 'R_Squared', 'Adj_R2', 'RMSE');
fprintf('------------------------------------------------------------------------\n');

for idx = 1:8
    if idx <= 4
        X_val = Z;
    else
        X_val = A;
    end
    Y_val = Epsilon_Matrix(:, idx);
    
    % Linear Fitting execution
    p = polyfit(X_val, Y_val, 1);
    y_fit = polyval(p, X_val);
    
    % Calculation of Residuals
    residuals = Y_val - y_fit;
    SS_res = sum(residuals.^2);
    SS_tot = sum((Y_val - mean(Y_val)).^2);
    
    % Statistical calculation formulas
    r_sq = 1 - (SS_res / SS_tot);
    n = length(X_val);
    adj_r_sq = 1 - ((1 - r_sq) * (n - 1) / (n - 2));
    rmse_val = sqrt(mean(residuals.^2));
    
    % Pearson Coefficient Extraction
    r_matrix = corrcoef(X_val, Y_val);
    r_coeff = r_matrix(1, 2);
    
    fprintf('Case %-5d %-12.4f %-10.4f %-10.4f %-10.4f\n', idx, r_coeff, r_sq, adj_r_sq, rmse_val);
end
fprintf('========================================================================\n');
%% 5. Define Global Uniform Scaling Bounds
global_y_min = min(Epsilon_Matrix(:)) - 5;
global_y_max = max(Epsilon_Matrix(:)) + 5;

z_min = min(Z) - 5; z_max = max(Z) + 5;
a_min = min(A) - 10; a_max = max(A) + 10;

%% 6. Plot Generation: Multi-Panel Comparative Matrix (4 for Z, 4 for A)
figure('Name', 'Heuristic Structural Optimization Matrix', 'Position', [50, 50, 1400, 850]);

titles = {'Case I: Standard (\epsilon vs Z)', 'Case II: \epsilon vs Z Variant', ...
          'Case III: \epsilon vs Z Variant B', 'Case IV: Inverted Num (\epsilon vs Z)', ...
          'Case V: \epsilon vs A Variant', 'Case VI: \epsilon vs A Variant B', ...
          'Case VII: Optimized Linear (\epsilon vs Atomic Mass)', 'Case VIII: Dual Sign Inversion (\epsilon vs A)'};

colors = {'#7E191B', '#D95319', '#EDB120', '#77AC30', '#4DBFFF', '#415A77', '#0072BD', '#A2142F'};

for idx = 1:8
    subplot(2, 4, idx);
    
    if idx <= 4
        plot(Z, Epsilon_Matrix(:, idx), 'o-', 'Color', colors{idx}, 'LineWidth', 1, 'MarkerSize', 3);
        xlim([z_min, z_max]);
        xlabel('Atomic Number (Z)', 'FontSize', 10, 'FontWeight', 'bold');
    else
        if idx == 7
            plot(A, Epsilon_Matrix(:, idx), 'o-', 'Color', colors{idx}, 'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', colors{idx});
        else
            plot(A, Epsilon_Matrix(:, idx), 'o-', 'Color', colors{idx}, 'LineWidth', 1, 'MarkerSize', 3);
        end
        xlim([a_min, a_max]);
        xlabel('Atomic Mass (A) [amu]', 'FontSize', 10, 'FontWeight', 'bold');
    end
    
    grid on;
    ylim([global_y_min, global_y_max]);
    ylabel('\epsilon (Error Parameter)', 'FontSize', 10, 'FontWeight', 'bold');
    title(titles{idx}, 'FontSize', 11, 'FontWeight', 'bold');
end

sgtitle('Heuristic Model Analysis: 4 Channels vs Z and 4 Channels vs Atomic Mass (Uniform Scaling)', 'FontSize', 16, 'FontWeight', 'bold');