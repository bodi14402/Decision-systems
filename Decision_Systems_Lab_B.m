clc; clear; close all;
%% Preliminaries
% In preparation for this laboratory, go to the ACS6124 Blackboard site and download
% the following files from the Decision Systems Laboratories area, and place them in
% your working directory (or other convenient location on the Matlab path):

% 1 • evaluateControlSystem.m and its utility functions getPeakInfo.m, getRiseTime.m,
% getSettlingTime.m, getShoots.m, and getSSError.m—these are the files needed
% to evaluate a control system design; PLACED IN "Evaluate_and_Utility" FILE

% 2 • sampling.zip—this is a set of open-source sampling plan routines and utilities 
% developed for engineering design problems by Alexander Forrester and colleagues 
% from the University of Southampton. (sampling file)

%% Add necessary directories to the path
% Get the current working directory
current_dir = pwd;
% Generate paths to the 'sampling' and 'Evaluate_and_Utility' directories
sampling_dir = genpath(fullfile(current_dir,'sampling'));
evaluate_utility_dir = genpath(fullfile(current_dir, 'Evaluate_and_Utility'));
% Add the directories and their subdirectories to the path
addpath(sampling_dir, evaluate_utility_dir);

%% Building the files provided
mex(fullfile(pwd, '/EA_Toolbox/rank_nds.c'));
mex(fullfile(pwd, '/EA_Toolbox/crowdingNSGA_II.c'));
mex(fullfile(pwd, '/EA_Toolbox/btwr.c'));
mex(fullfile(pwd, '/EA_Toolbox/sbx.c'));
mex(fullfile(pwd, '/EA_Toolbox/polymut.c'));
mex('-DVARIANT=4', ...
    fullfile(pwd, '/Hypervolume/Hypervolume_MEX.c'), ...
    fullfile(pwd, '/Hypervolume/hv.c'), ...
    fullfile(pwd, '/Hypervolume/avl.c'));
mex(fullfile(pwd, '/EA_Toolbox/rank_prf.c'));

% Add sampling and evaluation functions to the path
% Define the current directory
currentDir = pwd;

% Go up one level to the parent directory
parentDir = fileparts(currentDir);

% Add sampling and evaluation functions to the path
addpath(fullfile(parentDir, 'Lab_A', 'sampling'));
addpath(fullfile(parentDir, 'Lab_A', 'Evaluate_and_Utility'));

addpath(fullfile(pwd, '/EA_Toolbox/'));

%% Sampling Plans
%   _________                     .__  .__                 __________.__                        
%  /   _____/____    _____ ______ |  | |__| ____    ____   \______   \  | _____    ____   ______
%  \_____  \\__  \  /     \\____ \|  | |  |/    \  / ___\   |     ___/  | \__  \  /    \ /  ___/
%  /        \/ __ \|  Y Y  \  |_> >  |_|  |   |  \/ /_/  >  |    |   |  |__/ __ \|   |  \\___ \ 
% /_______  (____  /__|_|  /   __/|____/__|___|  /\___  /   |____|   |____(____  /___|  /____  >
%         \/     \/      \/|__|                \//_____/                       \/     \/     \/ 
% 

% Define design constraints and sampling plans using cell arrays
design_constraints = {'max pole', 'gain margin', 'phase margin', 'rise time', ...
                     'peak time', 'overshoot', 'undershoot', 'settling time', ...
                     'steady-state error', 'control input'};
sampling_plans_list = {'full factorial', 'sobol set', 'latin hypercube', 'random Latin hypercube'};

% Function handles for analysis and knowledge discovery
analyzeSamplingPlans = @(plans) mmphi_analysis(plans);
evaluateCS = @(P) evaluateControlSystem(P);
discoverKnowledge = @(Z, plan, constraints) knowledge_discovery(Z, plan, constraints);

% Perform analysis, evaluation, and knowledge discovery
[P, bestPlan] = analyzeSamplingPlans(sampling_plans_list);
fprintf('Best: %s\n', bestPlan);
Z = evaluateCS(P);
discoverKnowledge(Z, bestPlan, design_constraints);

PCA_and_Kmeans(evaluateControlSystem(P));
%% Building the optimizing engine

best_sampling_plan = [];
iterations = 2;
priority = [3, 2, 2, 1, 0, 1, 0, 0, 1, 2];
weighted_priority = [0.8, 0.6, 0.6, 0.4, 0.2, 0.4, 0.2, 0.2, 0.6, 0.6];

% Lab B goals
goals = [1, -6, 20, 2, 10, 10, 8, 20, 1, 0.67];
optimize_preferability(true, P, iterations, goals, priority, weighted_priority, best_sampling_plan, design_constraints);

% Assignment goals
% goals = [1, -6, 20, 2, 10, 10, 8, 20, 1, 0.63];
% buildOptimizingEngine(true, P, iterations, goals, priority, weighted_priority, best_sampling_plan, design_constraints);

%% Diffrent sampling plans input and outputs
    %                        .__    .__ 
    %   _____   _____ ______ |  |__ |__|
    %  /     \ /     \\____ \|  |  \|  |
    % |  Y Y  \  Y Y  \  |_> >   Y  \  |
    % |__|_|  /__|_|  /   __/|___|  /__|
    %       \/      \/|__|        \/    

% Inputs:
%       X - sampling plan
%       q - exponent used in the calculation of the metric
%       p - the distance metric to be used (p=1 rectangular - default, p=2
%           Euclidean)
%
% Output:
%       Phiq - sampling plan `space-fillingness' metric


    %        .__  .__     
    % _______|  | |  |__  
    % \_  __ \  | |  |  \ 
    %  |  | \/  |_|   Y  \
    %  |__|  |____/___|  /
    %                  \/ 

% Inputs:
%       n - desired number of points
%       k - number of design variables (dimensions)
%       Edges - if Edges=1 the extreme bins will have their centres on the
%               edges of the domain, otherwise the bins will be entirely 
%               contained within the domain (default setting). 
%
% Output:
%       X - Latin hypercube sampling plan of n points in k dimensions.


%   _____     .__  .__   _____               __               .__       .__    
% _/ ____\_ __|  | |  |_/ ____\____    _____/  |_  ___________|__|____  |  |   
% \   __\  |  \  | |  |\   __\\__  \ _/ ___\   __\/  _ \_  __ \  \__  \ |  |   
%  |  | |  |  /  |_|  |_|  |   / __ \\  \___|  | (  <_> )  | \/  |/ __ \|  |__ 
%  |__| |____/|____/____/__|  (____  /\___  >__|  \____/|__|  |__(____  /____/ 
%                                  \/     \/                          \/       

% Inputs:
%       q - k-vector containing the number of points along each dimension
%       Edges - if Edges=1 the points will be equally spaced from edge to
%               edge (default), otherwise they will be in the centres of 
%               n = q(1)*q(2)*...q(k) bins filling the unit cube.
%
% X - full factorial sampling plan

    %             ___.          .__                 __   
    %   __________\_ |__   ____ |  |   ______ _____/  |_ 
    %  /  ___/  _ \| __ \_/ __ \|  |  /  ___// __ \   __\
    %  \___ (  <_> ) \_\ \  ___/|  |__\___ \\  ___/|  |  
    % /____  >____/|___  /\___  >____/____  >\___  >__|  
    %      \/          \/     \/          \/     \/      
     
% p = sobolset(d) constructs a d-dimensional point set p, which is a 
% sobolset object with default property settings. The input argument d 
% corresponds to the Dimensions property of p.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
    % __________                   __  .__                      
    % \_  _____/_ __  ____   _____/  |_|__| ____   ____   ______
    %  |   __)|  |  \/    \_/ ___\   __\  |/  _ \ /    \ /  ___/
    %  |    \ |  |  /   |  \  \___|  | |  (  <_> )   |  \\___ \ 
    %  \___ / |____/|___|  /\___  >__| |__|\____/|___|  /____  >
    %     \/             \/     \/                    \/     \/ 

% Function to analyze different sampling plans using mmphi
% Inputs:
%   sampling_plans_list - a cell array of sampling plan names
% Outputs:
%   P_best - the best sampling plan based on mmphi metric
%   best_sampling_plan - name of the best sampling plan

function [P_best, best_sampling_plan] = mmphi_analysis(sampling_plans_list)
    % Initialize parameters
    params = struct('scale', 1, 'q', [10, 10], 'Edges', 1, 'phi_inf', Inf, ...
                'n_rows', ceil(length(sampling_plans_list) / 2), ...
                'best_sampling_plan', '', 'P_best', []);

    % Create a figure
    figure;

    % Loop over sampling plans
    for i = 1:length(sampling_plans_list)
        sampling_plan = strrep(sampling_plans_list{i}, ' ', '_');
        
        try
            % Generate the sampling plan
            switch sampling_plan
                case 'full_factorial'
                    X = fullfactorial(params.q, params.Edges);
                case 'sobol_set'
                    X = sobolset(length(params.q));
                    X = net(X, prod(params.q));
                case 'latin_hypercube'
                    X = lhsdesign(prod(params.q), length(params.q));
                case 'random_Latin_hypercube'
                    X = rlh(prod(params.q), length(params.q), params.Edges);
                otherwise
                    error('Invalid sampling plan specified.');
            end

            % Compute the MMPhi metric
            phi_metric = mmphi(X * params.scale, 5, 1);
            fprintf(' MMPhi for %s is: %f\n', sampling_plan, phi_metric);

            % Update the best sampling plan
            if phi_metric < params.phi_inf
                params.phi_inf = phi_metric;
                params.best_sampling_plan = sampling_plans_list{i};
                params.P_best = X;
            end

            % Plot the sampling plan
            plotSamplingPlan(X, params.n_rows, i, sampling_plans_list{i});
        catch ME
            fprintf('Error processing sampling plan %s: %s\n', sampling_plans_list{i}, ME.message);
        end
    end
    
    % Assign the best sampling plan and P_best for output
    best_sampling_plan = params.best_sampling_plan;
    P_best = params.P_best;
end

% % Function to update the best sampling plan
% function [min_phi_metric, best_sampling_plan, P_best] = update_best_sampling_plan(phi_metric, min_phi_metric, sampling_plan, X)
%     % Initialize outputs in case they are not assigned within the if statement
%     best_sampling_plan = '';
%     P_best = [];
% 
%     if phi_metric < min_phi_metric
%         min_phi_metric = phi_metric;
%         best_sampling_plan = sampling_plan;
%         P_best = X;
%     end
% end

%% Function to implement knowledge discovery
function knowledge_discovery(Z, best_sampling_plan, design_constraints)
    figure;

    % Create a parallel plot
    p = parallelplot(Z, 'Color', 'b');
    p.CoordinateTickLabels = design_constraints;
    p.YLabel = 'Performance metric value';
    p.Title = sprintf('Performance evaluations for %s sampling plan', best_sampling_plan);

    % Create a matrix scatter plot
    figure;
    plotmatrix(Z);
    title('Matrix Scatter Plot');

    % Create a grouped matrix scatter plot
    figure;
    gplotmatrix(Z);
    title('Grouped Matrix Scatter Plot');

    % Create a glyph plot
    figure;
    % glyphplot(Z, 'glyph','face','obslabels',design_constraints);
    glyphplot(Z)
    title('Glyph Plot');

    % % Perform principal component analysis
    % [coeff,score,~,~,~] = pca(Z);
    % figure;
    % biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',design_constraints);
    % title('Principal Component Analysis');

    % Perform k-means clustering
    idx = kmeans(Z,3);
    figure;
    gscatter(Z(:,1),Z(:,2),idx);
    title('K-Means Clustering');
end

%% Function to plot the sampling plan
function plotSamplingPlan(X, n_rows, i, sampling_plan)
    subplot(n_rows, 2, i);
    plot(X(:,1), X(:,2), 'bx'); 
    title(sprintf('%s', sampling_plan));
    xlabel('x_1');
    ylabel('x_2');
end

%% Function to optimize the sampling plan
function Z_optimized = optimizeControlSystem(P)
    % Evaluate the control system
    Z = evaluateControlSystem(P);

    % Optimization steps
    Z = convertGainMarginToDecibels(Z);
    Z = minimizeGainMargin(Z);
    Z = minimizeAbsoluteDeviation(Z);
    Z = removeInfValues(Z);
    Z_optimized = Z;
end

%% Helpers for Z_optimized
% Convert gain margin to decibels
function Z = convertGainMarginToDecibels(Z)
    Z(:,2) = 20*log10(Z(:,2));
end

% Minimize the gain margin
function Z = minimizeGainMargin(Z)
    Z(:,2) = -Z(:,2);
end

% Minimize the absolute deviation of the gain margin within 30 and 70 dB
function Z = minimizeAbsoluteDeviation(Z)
    Z(:,3) = abs(Z(:,3) - 50);
end

% Remove all inf values in Z, set to a high value
function Z = removeInfValues(Z)
    Z(isinf(Z)) = 1000;
end



%%


%% Function to analyze the best fit from the matrix of design evaluations, priority, and goals
function priority_based_analysis(decisionVars, evalMatrix, priorityLevels, weightedPriorities, targetGoals, constraintLabels)
    % Initialize priority level indices
    priorityIndices = cell(max(priorityLevels), 1);

    % Determine indices for each priority level
    for i = 1:max(priorityLevels)
        priorityIndices{i} = find(priorityLevels == i);
    end

    % Initialize level indices
    levelIndices = cell(max(priorityLevels), 1);

    % Get indices for each priority level
    for i = 1:max(priorityLevels)
        for j = 1:length(priorityIndices{i})
            levelIndices{i} = [levelIndices{i}; find(evalMatrix(:, priorityIndices{i}(j)) < targetGoals(priorityIndices{i}(j)))];
        end
    end

    % Get intersection of indices for all priority levels
    intersectIndices = levelIndices{1};
    for i = 2:max(priorityLevels)
        intersectIndices = intersect(intersectIndices, levelIndices{i});
    end

    % Determine final set of indices based on intersection results
    finalIndices = intersectIndices;
    if isempty(finalIndices)
        fprintf('No indices satisfy all constraints\n');
        for i = max(priorityLevels):-1:1
            if ~isempty(levelIndices{i})
                fprintf('The number of indices that satisfy the highest %d constraints is: %d\n', i, length(levelIndices{i}));
                finalIndices = levelIndices{i};
                break;
            end
        end
    else
        fprintf('The number of indices that satisfy all constraints is: %d\n', length(finalIndices));
    end

    % Check the difference between the target goals and the evaluation matrix values to find the best fit
    optimalFit = zeros(length(finalIndices), 1);
    for idx = 1:length(finalIndices)
        diff = abs(evalMatrix(finalIndices(idx), :) - targetGoals) ./ targetGoals;
        diff = diff .* weightedPriorities;
        % remove the NaN values
        diff(isnan(diff)) = 0;
        optimalFit(idx) = sum(diff);
    end

    % get the index of the minimum optimal fit
    minOptimalFitIdx = find(optimalFit == min(optimalFit));

    % print the row of evalMatrix that satisfies the minimum optimal fit
    optimalSolution = evalMatrix(finalIndices(minOptimalFitIdx), :);

   % print the violation of the constraints
    violatedConstraints = optimalSolution > targetGoals;
    if any(violatedConstraints)
        fprintf('The following constraints are NOT satisfied:\n');
        for i = find(violatedConstraints)
            fprintf('- %s\n', constraintLabels{i});
        end
    else
        fprintf('All constraints are satisfied.\n');
    end

    % print the optimal solution
    fprintf('The optimal solution is:\n');
    for i = 1:length(optimalSolution)
        fprintf('- %s: %f\n', constraintLabels{i}, optimalSolution(i));
    end

    % get the values of decision variables that satisfy the minimum optimal fit
    fprintf('The values of decision variables that satisfy the minimum optimal fit are:\n');
    for i = 1:size(decisionVars, 2)
        fprintf('- Variable %d: %f\n', i, decisionVars(finalIndices(minOptimalFitIdx), i));
    end
end

%% Function for optimizing with preferability
function optimize_preferability(enable_preference, P, iterations, goals, priority, weighted_priority, best_sampling_plan, design_constraints)
    % Initialize reference point and result
    reference_point = max(optimizeControlSystem(P));
    Res = [];
    % figure;

    for i = 1:iterations

        Z = initialize_population(P);
        [ranking, distances] = calculate_fitness(Z, enable_preference, goals, priority);
        selectThese = perform_selection(distances);
        unifiedPop = perform_variation(P, selectThese);
        P = perform_survival(unifiedPop, Z, enable_preference, goals, priority);
        Res = check_convergence(Z, reference_point, Res);
        plotProgress(P, i, Res)

    end

    Z = optimizeControlSystem(P);
    knowledge_discovery(Z, best_sampling_plan, design_constraints);
    priority_based_analysis(P, Z, priority, weighted_priority, goals, design_constraints);
end

%% Helper
% Initialize population
function Z = initialize_population(P)
    Z = optimizeControlSystem(P);
end

% Calculate fitness
function [ranking, distances] = calculate_fitness(Z, enable_preference, goals, priority)
    if enable_preference
        ranking = rank_prf(Z, goals, priority);
    else
        ranking = rank_nds(Z);
    end
    ranking = max(ranking) - ranking;
    distances = crowding(Z, ranking);
end

% Perform selection
function selectThese = perform_selection(distances)
    selectThese = btwr(distances, length(distances));
end

% Perform variation
function unifiedPop = perform_variation(P, selectThese)
    parents  = P(selectThese, :);
    bounds = [0, 0; 1, 1];
    offspring = sbx(parents, bounds);
    C = polymut(offspring, bounds);
    unifiedPop = [P; C];
end

% Perform survival
function P = perform_survival(unifiedPop, ~, enable_preference, goals, priority)
    Z_unified = optimizeControlSystem(unifiedPop);
    if enable_preference
        new_indices = reducerNSGA_II(unifiedPop, rank_prf(Z_unified, goals, priority), crowding(Z_unified, rank_prf(Z_unified, goals, priority)));
    else
        new_indices = reducerNSGA_II(unifiedPop, rank_nds(Z_unified), crowding(Z_unified, rank_nds(Z_unified)));
    end
    P = unifiedPop(new_indices, :);
end


% PCA for Appropriate data mining and visualisation 
function PCA_and_Kmeans(Z)
    % Preprocess the data and apply PCA
    score = preprocess_and_pca(Z);

    % Perform kmeans clustering and plot the results
    kmeans_and_plot(score);
end

% Function to preprocess the data and apply PCA
function score = preprocess_and_pca(Z)
    % Remove inf values and normalize the data
    Z(isinf(Z)) = 1e6;
    Z = normalize(Z);

    % Apply PCA and return the score
    [~, score] = pca(Z);
end

% Function to perform kmeans clustering and plot the results
function kmeans_and_plot(score)
    % Perform kmeans clustering on the first three principal components
    [idx, C] = kmeans(score(:,1:3), 3);

    % Plot the first three principal components in a new figure
    figure;
    scatter3(score(:,1), score(:,2), score(:,3), 50, idx, 'filled', 'MarkerEdgeColor', 'k');
    hold on;
    scatter3(C(:,1), C(:,2), C(:,3), 200, 'r', 'filled', 'MarkerEdgeColor', 'k');
    hold off;
    title('Principal Components 1-3');
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    legend('Clusters', 'Centroids');
    grid on;

    % Perform kmeans clustering on the next three principal components
    [idx, C] = kmeans(score(:,4:6), 3);

    % Plot the next three principal components in a new figure
    figure;
    scatter3(score(:,4), score(:,5), score(:,6), 50, idx, 'filled', 'MarkerEdgeColor', 'k');
    hold on;
    scatter3(C(:,1), C(:,2), C(:,3), 200, 'r', 'filled', 'MarkerEdgeColor', 'k');
    hold off;
    title('Principal Components 4-6');
    xlabel('PC4');
    ylabel('PC5');
    zlabel('PC6');
    legend('Clusters', 'Centroids');
    grid on;
    pause(0.5);
end

% Check for convergence
function output_results = check_convergence(Z, reference_point, output_results)
    output_results = [output_results, Hypervolume_MEX(Z, reference_point)];
end

% Plot progress
function plotProgress(P, i, Res, ~)
        % clf(figure(15));
        figure(15);
        plot(P(:,1), P(:,2), 'x');
        title(sprintf('Sampling Plan Update %d', i));
        xlabel('x_1');
        ylabel('x_2');
        drawnow;
        
        % clf(figure(16));
        figure(16);
        plot(Res, 'LineWidth', 2, 'Color', 'b');
        title('Hypervolume Convergence');
        xlabel('Iteration');
        ylabel('Hypervolume');
        xlim([0, 51]);
        ylim([1e55, 0.5* 1e56]);
        drawnow;
end

