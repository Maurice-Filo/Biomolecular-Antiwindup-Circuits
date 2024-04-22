function [T, X] = SSA_Grid(Propensity, Stoichiometry, NetworkParameters, X0, TimeSpan, N_Grid)
% SSA simulates a chemical reaction network using the Direct Method of the Gillespie Stochastic Simulation Algorithm.
% The simulation is stored on a predefined grid.
%
% Inputs:
% - Propensity: Function handle to calculate propensities.
% - Stoichiometry: Stoichiometric coefficients of the reactions.
% - NetworkParameters: Parameters needed for the propensity function.
% - X0: Initial state.
% - TimeSpan: Vector [t0, tf] where t0 is the initial time and tf is the final time.
% - N_Grid: Number of grid points.
%
% Outputs:
% - T: Regularly spaced time grid.
% - X: State of the system on the time grid T.
t0 = TimeSpan(1);
tf = TimeSpan(2);
N_Species = length(X0);
T = linspace(t0, tf, N_Grid);
dT = T(2) - T(1);
X = [X0, zeros(N_Species, N_Grid-1)];
X_Current = X0;
T_Current = T(1);
Index = 1;
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
    T_Current = T_Current + tau;
    X_Current = X_Current + Stoichiometry(:,j);
    % Store
    if T_Current >= T(Index + 1)
        Index_New = floor(T_Current/dT + 1);
        X(:,Index+1:min(Index_New, N_Grid)) = repmat(X_Current, 1, min(Index_New - Index, N_Grid - Index));
        Index = Index_New;
    end
end
end

