function [T, X] = SSA(Propensity, Stoichiometry, NetworkParameters, X0, TimeSpan, N_MaxRxnEvents)
% SSA simulates a chemical reaction network using the Direct Method of the Gillespie Stochastic Simulation Algorithm.
% The sizes of the time and state vectors are preallocated.
%
% Inputs:
% - Propensity: Function handle to calculate propensities.
% - Stoichiometry: Stoichiometric coefficients of the reactions.
% - NetworkParameters: Parameters needed for the propensity function.
% - X0: Initial state.
% - TimeSpan: Vector [t0, tf] where t0 is the initial time and tf is the final time.
% - N_MaxRxnEvents: Maximum number of reaction events to simulate.
%
% Outputs:
% - T: Time points at which reactions occur.
% - X: State of the system at the time points T.

t0 = TimeSpan(1);
tf = TimeSpan(2);
N_Species = length(X0);
T = zeros(1, N_MaxRxnEvents);
X = [X0, zeros(N_Species, N_MaxRxnEvents)];
RxnCounter = 1;
T(1) = t0;
while T(RxnCounter) < tf
    % Calculate propensity at the current (x(t), t)
    a = Propensity(X(:,RxnCounter), NetworkParameters); 
    % Generate a random time at which the next reaction will occur
    a0 = sum(a);
    r1 = rand(1);
    tau = -log(r1) / a0;  
    % Randomly identify which reaction will occur
    r2 = rand(1);
    j = find((cumsum(a) >= r2*a0), 1);
    % Carry out the reaction and update time
    RxnCounter = RxnCounter + 1;
    T(RxnCounter) = T(RxnCounter - 1) + tau;
    X(:, RxnCounter) = X(:, RxnCounter - 1) + Stoichiometry(:,j);
    if RxnCounter >= N_MaxRxnEvents
        break
    end        
end
T(RxnCounter:end) = [];
X(:, RxnCounter:end) = [];
end

