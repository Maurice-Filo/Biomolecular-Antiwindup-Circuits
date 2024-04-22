function [t_vector, MeanX, VarX] = ComputeStatistics(T, X, N)
% ComputeStatistics computes the mean and variance trajectories of multiple species over 
% a shared discretized time horizon. The function resamples the given trajectories using 
% a zero-order hold method between jumping times.
%
% Inputs:
% - T: Cell array containing time vectors for multiple trajectories.
% - X: Cell array containing state matrices for multiple trajectories.
% - N: Number of samples for the discretized time horizon.
%
% Outputs:
% - t_vector: Discretized time horizon.
% - MeanX: Matrix of mean trajectories for all species.
% - VarX: Matrix of variance trajectories for all species.

    %% Find the maximum time among the different stochastic trajectories
    t_max = 1e20;
    for i = 1 : length(T)
        T_Vector = T{i};
        t_max = min(t_max, T_Vector(end));
    end
    
    %% Discretize the Time Horizon
    t_vector = linspace(0, t_max, N);
        
    %% Compute the sample paths at the common time horizon
    X_Resampled = cell(length(T), 1);
    for i = 1 : length(T)
        [X_Resampled{i}] = ResampleZOH(X{i}, T{i}, t_vector);
    end
    
    %% Compute the Statistics
    MeanX = zeros(size(X{1}, 1), length(t_vector));
    VarX = zeros(size(X{1}, 1), length(t_vector));
    for j = 1 : size(X{1}, 1)
        Xj = zeros(length(T), length(t_vector));
        for i = 1 : length(T)
            Xj(i,:) = X_Resampled{i}(j,:);
        end
        MeanX(j,:) = mean(Xj, 1);
        VarX(j,:) = var(Xj, 1);
    end
end

