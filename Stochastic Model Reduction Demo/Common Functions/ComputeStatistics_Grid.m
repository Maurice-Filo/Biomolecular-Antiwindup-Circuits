function [MeanX, VarX] = ComputeStatistics_Grid(T, X)
% ComputeStatistics_Grid computes the mean and variance trajectories of multiple 
% species over a uniformly gridded time horizon.
%
% Inputs:
% - T: Cell array containing uniformly gridded time vectors for multiple trajectories.
% - X: Cell array containing state matrices for multiple trajectories.
%
% Outputs:
% - MeanX: Matrix of mean trajectories for all species.
% - VarX: Matrix of variance trajectories for all species.

    %% Compute the Statistics
    MeanX = zeros(size(X{1}, 1), length(T{1}));
    VarX = zeros(size(X{1}, 1), length(T{1}));
    for j = 1 : size(X{1}, 1)
        Xj = zeros(length(T), length(T{1}));
        for i = 1 : length(T)
            Xj(i,:) = X{i}(j,:);
        end
        MeanX(j,:) = mean(Xj, 1);
        VarX(j,:) = var(Xj, 1);
    end
end

