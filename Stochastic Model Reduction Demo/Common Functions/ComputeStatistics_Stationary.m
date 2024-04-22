function [StationaryMean, StationaryVariance, SimTime] = ComputeStatistics_Stationary(PropensityFunction, StoichiometryMatrix, Parameters, X0, TimeSpan, TransientPortion, N_Trajectories)
% ComputeStatistics_Stationary estimates the stationary mean and variance
% of a stochastic system.
%
% Inputs:
% - Propensity: A function handle to compute propensities.
% - Stoichiometry: Stoichiometry matrix of the reactions.
% - NetworkParameters: Parameters specific to the network.
% - X0: Initial state of the system.
% - TimeSpan: A 2-element vector specifying the start and end time.
% - TransientPortion: Portion of the simulation treated as transient and discarded for statistics.
% - N_Trajectories: Number of trajectories/samples for statistics computation.
%
% Outputs:
% - StationaryMean: Estimated stationry mean.
% - StationaryVariance: Estimated stationary variance.
% - SimTime: Total execution time of the simulation.

%% Initialize simulation timer
TStart = tic;

%% Preallocate arrays for efficiency
MeanX = zeros(length(X0), N_Trajectories);
VarX = zeros(length(X0), N_Trajectories);

%% Compute statistics using multiple trajectories in parallel
parfor i = 1 : N_Trajectories
    [MeanX(:,i), VarX(:,i)] = SSA_Stationary(PropensityFunction, StoichiometryMatrix, Parameters, X0, TimeSpan, TransientPortion);
end

%% Calculate the mean and variance across all trajectories
StationaryMean = mean(MeanX, 2);
StationaryVariance = mean(VarX, 2);

%% Get simulation time
SimTime = toc(TStart);
end

