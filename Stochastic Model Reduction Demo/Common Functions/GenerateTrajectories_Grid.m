function [T, X, SimTime] = GenerateTrajectories_Grid(StoichiometryMatrix, PropensityFunction, Parameters, IC, TimeSpan, N_Trajectories, N_Grid)
% GenerateTrajectories_Grid generates multiple trajectories for a chemical reaction system using grid-based SSA.
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
% - T: Cell array where each entry is a time vector for a trajectory.
% - X: Cell array where each entry is a state matrix for a trajectory.
% - SimTime: Total time taken for the simulations.

TStart = tic;
T = cell(N_Trajectories,1);
X = cell(N_Trajectories,1);
parfor i = 1 : N_Trajectories
    [T{i}, X{i}] = SSA_Grid(PropensityFunction, StoichiometryMatrix, Parameters, IC, TimeSpan, N_Grid);
end
SimTime = toc(TStart);
end