function [MeanX, VarX] = SSA_Stationary(Propensity, Stoichiometry, NetworkParameters, X0, TimeSpan, TransientPortion)
% SSA_Stationary estimates the stationary mean and variance of a stochastic
% system using a single trajectory.
% 
% Inputs:
% - Propensity: A function handle to compute propensities
% - Stoichiometry: Stoichiometry matrix of the reactions
% - NetworkParameters: Parameters of the network
% - X0: Initial state of the system
% - TimeSpan: A 2-element vector specifying the start and end time
% - TransientPortion: Portion of the simulation that should be treated as transient and discarded for statistics
%
% Outputs:
% - MeanX: Estimated stationary mean
% - VarX: Estimated stationary variance


%% Define time boundaries
t0 = TimeSpan(1);
tf = TimeSpan(2);
ts = t0 + TransientPortion * (tf - t0);

%% Initialize state variables
X_Current = X0;
T_Current = t0;
MeanX = 0*X0;
MeanX2 = 0*X0;

%% Simulate
while T_Current < tf
    % Calculate propensity at the current (x(t), t)
    a = Propensity(X_Current, NetworkParameters); 
    % Generate a random time at which the next reaction will occur
    a0 = sum(a);
    r1 = rand(1);
    tau = -log(r1) / a0;  
    % Randomly identify which reaction will occur
    r2 = rand(1);
    j = find((cumsum(a) >= r2*a0), 1);
    % Carry out the reaction and update time
    X_Old = X_Current;
    T_Current = T_Current + tau;
    X_Current = X_Old + Stoichiometry(:,j);
    % Compute moments if past the transient time
    if T_Current > ts
        MeanX = MeanX + (X_Old + X_Current) * tau / 2;
        MeanX2 = MeanX2 + (X_Old.^2 + X_Current.^2) * tau / 2;
    end    
end

%% Compute final statistics
MeanX = MeanX / (tf - ts);
VarX = MeanX2 / (tf - ts) - MeanX.^2;
end

