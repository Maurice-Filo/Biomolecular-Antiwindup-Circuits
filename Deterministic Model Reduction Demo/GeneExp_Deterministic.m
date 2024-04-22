function [t, X, X_Reduced_Mapped, Setpoint, Error] = GeneExp_Deterministic(eta, gamma_1, gamma_2, tf, t_Disturbance, DisturbanceFactor)
%% Controller 
% Controller Selection
StoichiometryMatrix_Controller = StoichiometryMatrix_AIF();
PropensityFunction_Controller = @PropensityFunction_AIF;
% Parameters
Parameters_Controller.mu = 10;
Parameters_Controller.eta = eta;
Parameters_Controller.delta_1 = 0;
Parameters_Controller.delta_2 = 0;
M = 2;                  % Number of Controller Species
% Actuation
Parameters_Controller.ActuationParameters.k = 1;
Parameters_Controller.ActuationFunction = @(z_1, Parameters) (Parameters.k * z_1);
% Sensing
Parameters_Controller.SensingParameters.theta = 1;
Parameters_Controller.SensingFunction = @(y, Parameters) (Parameters.theta * y);

%% Reduced Controller
% Controller Selection
StoichiometryMatrix_Controller_Reduced = StoichiometryMatrix_AIF_Reduced();
PropensityFunction_Controller_Reduced = @PropensityFunction_AIF_Reduced; 

%% Plant
% Plant Selection
StoichiometryMatrix_Plant = StoichiometryMatrix_GeneExp();
PropensityFunction_Plant = @PropensityFunction_GeneExp;
% Parameters
Parameters_Plant.k_0 = 0;
Parameters_Plant.k_1 = 1;
Parameters_Plant.gamma_1 = gamma_1;
Parameters_Plant.gamma_2 = gamma_2;
L = 2;                  % Number of Plant Species
Input_Index = 1;
Output_Index = 2;

%% Simulation Settings
N_t = 1000;             % Time Samples
Solver = 'ODE15s';      % ODE Solver
% Initial Conditions
IC = zeros(L+M,1);      
IC_Reduced = [IC(1:L); IC(L+M) - IC(L+M-1)];    

%% Disturbance 
DisturbanceFactor_1 = DisturbanceFactor(1);
t_Disturbance_1 = t_Disturbance(1);
DisturbanceFactor_2 = DisturbanceFactor(2);
t_Disturbance_2 = t_Disturbance(2);
DisturbedParameter = 'mu';

%% Closed-Loop Network
% Parameters                    
Parameters_CL.Controller = Parameters_Controller;
Parameters_CL.Plant = Parameters_Plant;                  

% Stoichiometry Matrix
S_Mutual = zeros(L, size(StoichiometryMatrix_Controller,2));
S_Mutual(Input_Index, 1) = 1;
StoichiometryMatrix_CL = [	StoichiometryMatrix_Plant,                          S_Mutual; ...
                            zeros(M, size(StoichiometryMatrix_Plant,2)),        StoichiometryMatrix_Controller; ...
                         ];
% Propensity Function
PropensityFunction_CL = @(X, Parameters_CL) ...
                    ([ PropensityFunction_Plant(X(1:L), Parameters_CL.Plant); ...
                       PropensityFunction_Controller(X(L+1:L+M), X(L), Parameters_CL.Controller); ...
                    ]);
                
%% Reduced Closed-Loop Network
% Stoichiometry Matrix
S_Mutual = zeros(L, size(StoichiometryMatrix_Controller_Reduced,2));
S_Mutual(Input_Index, 1) = 1;
StoichiometryMatrix_CL_Reduced = [	StoichiometryMatrix_Plant,                                  S_Mutual; ...
                                    zeros(1, size(StoichiometryMatrix_Plant,2)),                StoichiometryMatrix_Controller_Reduced; ...
                                 ];
% Propensity Function
PropensityFunction_CL_Reduced = @(X, Parameters_CL) ...
                                ([ PropensityFunction_Plant(X(1:L), Parameters_CL.Plant); ...
                                PropensityFunction_Controller_Reduced(X(L+1:L+1), X(L), Parameters_CL.Controller); ...
                                ]);
                
%% Simulating the Full & Reduced Models
Setpoint = Parameters_CL.Controller.mu / Parameters_CL.Controller.SensingParameters.theta;
[t_1, X_1] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, IC, t_Disturbance_1, N_t, Solver);
[~, X_1_Reduced] = DSA(StoichiometryMatrix_CL_Reduced, PropensityFunction_CL_Reduced, Parameters_CL, IC_Reduced, t_Disturbance_1, N_t, Solver);
Parameters_CL.Controller.(DisturbedParameter) = Parameters_CL.Controller.(DisturbedParameter) * DisturbanceFactor_1;
Setpoint = [Setpoint, Parameters_CL.Controller.mu / Parameters_CL.Controller.SensingParameters.theta];
[t_2, X_2] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, X_1(:,end), t_Disturbance_2 - t_Disturbance_1, N_t, Solver);
[~, X_2_Reduced] = DSA(StoichiometryMatrix_CL_Reduced, PropensityFunction_CL_Reduced, Parameters_CL, X_1_Reduced(:,end), t_Disturbance_2 - t_Disturbance_1, N_t, Solver);
Parameters_CL.Controller.(DisturbedParameter) = Parameters_CL.Controller.(DisturbedParameter) * DisturbanceFactor_2;
Setpoint = [Setpoint, Parameters_CL.Controller.mu / Parameters_CL.Controller.SensingParameters.theta];
[t_3, X_3] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, X_2(:,end), tf - t_Disturbance_2, N_t, Solver);
[~, X_3_Reduced] = DSA(StoichiometryMatrix_CL_Reduced, PropensityFunction_CL_Reduced, Parameters_CL, X_2_Reduced(:,end), tf - t_Disturbance_2, N_t, Solver);
t = [t_1, t_Disturbance_1 + t_2(2:end), t_Disturbance_2 + t_3(2:end)];
X = [X_1, X_2(:,2:end), X_3(:,2:end)];
X_Reduced = [X_1_Reduced, X_2_Reduced(:,2:end), X_3_Reduced(:,2:end)];
X_Reduced_Mapped = [X_Reduced(1:L, :); max(X_Reduced(L+1, :), 0); max(-X_Reduced(L+1, :), 0)];
Error = norm(X - X_Reduced_Mapped, 'fro') / norm(X, 'fro');
end

