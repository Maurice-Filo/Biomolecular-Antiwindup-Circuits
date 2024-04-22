%% Clear Workspace
close all
clear
clc

%% Controller 
% Controller Selection
StoichiometryMatrix_Controller = StoichiometryMatrix_AIF();
PropensityFunction_Controller = @PropensityFunction_AIF;
% Parameters
Parameters_Controller.mu = 10;
Parameters_Controller.eta = 1e5;
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
Parameters_Plant.gamma_1 = 0.5;
Parameters_Plant.gamma_2 = 0.5;
L = 2;                  % Number of Plant Species
Input_Index = 1;
Output_Index = 2;

%% Simulation Settings
tf = 300;               % Final Time
N_t = 1000;             % Time Samples
Solver = 'ODE15s';      % ODE Solver
% Initial Conditions
IC = zeros(L+M,1);      
IC_Reduced = [IC(1:L); IC(L+M) - IC(L+M-1)];    

%% Disturbance 
DisturbanceFactor_1 = 2;
t_Disturbance_1 = 100;
DisturbanceFactor_2 = 1/3;
t_Disturbance_2 = 150;
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
Start = tic;
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
Error = norm(X - X_Reduced_Mapped, 'fro') / norm(X, 'fro')
SimTime = toc(Start)

%% Figure Settings
ScalingFactor = 1;
SS = 4;
Figure_Width = 12 * SS;
Figure_Height = 12 * SS;
FontSize = ScalingFactor*5 * SS;
FontSize_Small = ScalingFactor*3 * SS;
FontSize_Large = ScalingFactor*6 * SS;
LineWidth = ScalingFactor*0.65 * SS;
LineWidth_Thick = ScalingFactor * SS;
LineWidth_Thin = ScalingFactor*0.1 * SS;
MarkerSize = ScalingFactor*2 * SS;
Opacity = 0.5;
Colors = [39, 52, 139; ... % Blue
          227, 6, 19; ... % Red
          0, 152, 58; ... % Green
          230, 0, 126; ... % Magenta
          ]/255;
 
%% Set Figure 1
Handle_Figure1 = figure();
    Handle_Figure1.Color = [1 1 1];
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [0, 0, Figure_Width, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];

%% Axis for Simulations
Plot_Names = {'x_1', 'x_2', 'z_1', 'z_2'};
rows = 2;
columns = 2;
Handle_Axis1 = gobjects(rows, columns);
for i = 1 : rows
    for j = 1 : columns
        Handle_Axis1(i,j) = subplot(rows, columns, (i - 1)*columns + j); 
        Handle_Axis1(i,j).Box = 'on';
        Handle_Axis1(i,j).BoxStyle = 'full';
        Handle_Axis1(i,j).LineWidth = LineWidth_Thin;
        Handle_Axis1(i,j).FontSize = FontSize;
        hold(Handle_Axis1(i,j), 'on');
        grid(Handle_Axis1(i,j), 'on');
        Handle_Axis1(i,j).XMinorGrid = 'on';
        Handle_Axis1(i,j).YMinorGrid = 'on';
        Handle_Axis1(i,j).XLabel.String = 't';
        Handle_Axis1(i,j).YLabel.String = Plot_Names{(i - 1)*columns + j};
        plot(Handle_Axis1(i,j), t, X((i - 1)*columns + j,:), 'Color', [Colors(1,:), Opacity], 'LineWidth', 1.5*LineWidth_Thick, 'LineStyle', '-');
        plot(Handle_Axis1(i,j), t, X_Reduced_Mapped((i - 1)*columns + j,:), 'Color', Colors(2,:), 'LineWidth', LineWidth_Thick/3, 'LineStyle', '-');
        if (i - 1)*columns + j == L
            plot(Handle_Axis1(i,j), [0, t_Disturbance_1, t_Disturbance_1, t_Disturbance_2, t_Disturbance_2, tf], ...
                                    [Setpoint(1), Setpoint(1), Setpoint(2), Setpoint(2), Setpoint(3), Setpoint(3)], ...
                                    'LineWidth', LineWidth, 'Color', 'k', 'LineStyle', '--');
        end
    end
end
    
%% Save Figures
Save_Flag = 0;
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
    print(Handle_Figure1, 'NE_Response', '-dpdf', '-painters');
    Handle_Figure1.Color = [1, 1, 1];
end
