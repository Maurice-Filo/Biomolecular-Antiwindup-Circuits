function [T, Mean, Covariance, Variance, Skewness, SimTime] = GenerateMoments(StoichiometryMatrix, PropensityFunction, Parameters, IC, TimeSpan, N_Trajectories, N_Grid)
% GenerateMoments generates the dynamics of the first few moments based on multiple trajectories for a chemical reaction system using grid-based SSA.
%
% Inputs:
% - StoichiometryMatrix: Stoichiometric coefficients of the reactions.
% - PropensityFunction: Function handle to calculate propensities.
% - Parameters: Parameters needed for the propensity function.
% - IC: Initial state of the system.
% - TimeSpan: Vector [t0, tf] where t0 is the initial time and tf is the final time.
% - N_Trajectories: Number of trajectories to generate.
% - N_Grid: Number of grid points in the time domain.
%
% Outputs:
% - T: Time Grid
% - M: Array where each row is a moment trajectory.
% - SimTime: Total time taken for the simulations.

TStart = tic;
N_Species = size(StoichiometryMatrix, 1);
T = linspace(TimeSpan(1), TimeSpan(2), N_Grid);
M_1 = zeros(N_Species, N_Grid);
M_2 = zeros(N_Species, N_Species, N_Grid);
M_3 = zeros(N_Species, N_Grid);
parfor i = 1 : N_Trajectories
    [~, X] = SSA_Grid(PropensityFunction, StoichiometryMatrix, Parameters, IC, TimeSpan, N_Grid); 
    M_1 = M_1 + X;
    M_2 = M_2 + bsxfun(@times, reshape(X, [N_Species, 1, N_Grid]), reshape(X, [1, N_Species, N_Grid]));
    M_3 = M_3 + X.^3;
end
M_1 = M_1 / N_Trajectories;
M_2 = M_2 / N_Trajectories;
M_3 = M_3 / N_Trajectories;
Mean = M_1;
Covariance = M_2 - bsxfun(@times, reshape(Mean, [N_Species, 1, N_Grid]), reshape(Mean, [1, N_Species, N_Grid]));
Variance = reshape(Covariance, [N_Species^2, N_Grid]);
Variance = Variance(1:N_Species+1:end,:);
Skewness = (M_3 - 3* M_1 .* Variance - M_1.^3) ./ sqrt(Variance).^3;
SimTime = toc(TStart);
end