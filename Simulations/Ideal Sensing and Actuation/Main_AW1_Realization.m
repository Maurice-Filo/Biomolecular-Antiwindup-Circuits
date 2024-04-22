%% Clear Workspace
close 
clear
clc

%% Controller 
% Controller Selection
StoichiometryMatrix_Controller = StoichiometryMatrix_AW1_Realization();
PropensityFunction_Controller = @PropensityFunction_AW1_Realization;
% Parameters
Parameters_Controller.mu = 10;
Parameters_Controller.eta = 100;
Parameters_Controller.gamma_1 = 0;
Parameters_Controller.gamma_2 = 0;
M = 6;                  % Number of Controller Species
% Actuation
Parameters_Controller.ActuationParameters.k = 1;
Parameters_Controller.ActuationFunction = @(z_1, Parameters) (Parameters.k * z_1);
% Sensing
Parameters_Controller.SensingParameters.theta = 1;
Parameters_Controller.SensingFunction = @(y, Parameters) (Parameters.theta * y);
% Anti-Windup
Parameters_Controller.eta_P = 100;
Parameters_Controller.eta_PP = 100;
Parameters_Controller.beta_P = 5;
Parameters_Controller.beta_PP = 5;
Parameters_Controller.alpha_1 = 1;
Parameters_Controller.alpha_2 = 1;
Parameters_Controller.kappa_1 = 10;
Parameters_Controller.kappa_2 = 10;
Parameters_Controller.gamma_1_P = 1;
Parameters_Controller.gamma_2_P = 1;
Parameters_Controller.gamma_1_PP = 1;
Parameters_Controller.gamma_2_PP = 1;
Parameters_Controller.g_1 = @(z_1_PP, z_2_PP) Parameters_Controller.alpha_1*z_2_PP / (1 + z_1_PP/Parameters_Controller.kappa_1);
Parameters_Controller.g_2 = @(z_1_P, z_2_P) Parameters_Controller.alpha_2*z_2_P / (1 + z_1_P/Parameters_Controller.kappa_2);
Parameters_Controller.h_1 = @(z_1) z_1;
Parameters_Controller.h_2 = @(z_2) z_2;

%% Plant
% Plant Selection
StoichiometryMatrix_Plant = StoichiometryMatrix_GeneExp();
PropensityFunction_Plant = @PropensityFunction_GeneExp;
% Parameters
Parameters_Plant.k_0 = 6;
Parameters_Plant.k_1 = 1;
Parameters_Plant.gamma_1 = 1;
Parameters_Plant.gamma_2 = 1;
L = 2;                  % Number of Plant Species
Input_Index = 1;
Output_Index = 2;

%% Simulation Settings
tf = 250;               % Final Time
N_t = 1000;             % Time Samples
Solver = 'ODE23s';      % ODE Solver
IC = zeros(L+M,1);      % Initial Condition

%% Disturbance 
DisturbanceFactor_1 = 2;
t_Disturbance_1 = 100;
DisturbanceFactor_2 = 0.5;
t_Disturbance_2 = 150;
DisturbedParameter = 'k_0';

%% Closed-Loop Network
% Parameters                    
Parameters_CL.Controller = Parameters_Controller;
Parameters_CL.Plant = Parameters_Plant;                  

% Stoichiometry Matrix
S_Mutual = zeros(L, size(StoichiometryMatrix_Controller,2));
S_Mutual(1,Input_Index) = 1;
StoichiometryMatrix_CL = [	StoichiometryMatrix_Plant,                          S_Mutual; ...
                            zeros(M, size(StoichiometryMatrix_Plant,2)),        StoichiometryMatrix_Controller; ...
                         ];
% Propensity Function
PropensityFunction_CL = @(X, Parameters_CL) ...
                    ([ PropensityFunction_Plant(X(1:L), Parameters_CL.Plant); ...
                       PropensityFunction_Controller(X(L+1:L+M), X(L), Parameters_CL.Controller); ...
                    ]);
                
%% Simulating the Full Model
[t_1, X_1] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, IC, t_Disturbance_1, N_t, Solver);
Parameters_CL.Plant.(DisturbedParameter) = Parameters_CL.Plant.(DisturbedParameter) * DisturbanceFactor_1;
[t_2, X_2] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, X_1(:,end), t_Disturbance_2 - t_Disturbance_1, N_t, Solver);
Parameters_CL.Plant.(DisturbedParameter) = Parameters_CL.Plant.(DisturbedParameter) * DisturbanceFactor_2;
[t_3, X_3] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, X_2(:,end), tf - t_Disturbance_2, N_t, Solver);
t = [t_1, t_Disturbance_1 + t_2(2:end), t_Disturbance_2 + t_3(2:end)];
X = [X_1, X_2(:,2:end), X_3(:,2:end)];

%% Figure Settings
ScalingFactor = 1;
SS = 4;
Figure_Width = 6.5 * SS;
Figure_Height = 5 * SS;
FontSize = ScalingFactor*5 * SS;
FontSize_Small = ScalingFactor*3 * SS;
FontSize_Large = ScalingFactor*6 * SS;
LineWidth = ScalingFactor*0.65 * SS;
LineWidth_Thick = ScalingFactor*0.8 * SS;
LineWidth_Thin = ScalingFactor*0.1 * SS;
MarkerSize = ScalingFactor*5 * SS;
Opacity = 0.5;
Colors = lines(10);
    
%% Set Figure 1
Handle_Figure1 = figure();
    Handle_Figure1.Color = [1 1 1];
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [0, 0, Figure_Width, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];
    
%% Set Figure 2
Handle_Figure2 = figure();
    Handle_Figure2.Color = [1 1 1];
    Handle_Figure2.PaperUnits = 'centimeters';
    Handle_Figure2.Units = 'centimeters';
    Handle_Figure2.Position = [0, 0, Figure_Width, Figure_Height];
    Handle_Figure2.PaperPositionMode = 'auto';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)];

%% Axis for Simulations
Plot_Indeces = [L, L+1, L+2];
Plot_Names = {'$x_L$', '$z_1$', '$z_2$'};

Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Box = 'on';
  	Handle_Axis1.BoxStyle = 'full';
 	Handle_Axis1.LineWidth = LineWidth_Thin;
  	Handle_Axis1.FontSize = FontSize;
 	hold(Handle_Axis1, 'on');
 	grid(Handle_Axis1, 'on');
   	Handle_Axis1.XMinorGrid = 'on';
  	Handle_Axis1.YMinorGrid = 'on';
   	Handle_Axis1.Title.String = 'Integral Windup';
  	Handle_Axis1.XTickLabel = [0, t_Disturbance_1, t_Disturbance_2, tf];
  	Handle_Axis1.XLabel.String = 'Time';
    Handle_Axis1.YLabel.String = 'Response';
for i = 1 : length(Plot_Indeces)
    plot(Handle_Axis1, t, X(Plot_Indeces(i), :), 'Color', Colors(i,:), 'LineWidth', LineWidth_Thick);
end
Handle_Legend = legend(Handle_Axis1, Plot_Names);
    Handle_Legend.Interpreter = 'latex';
    Handle_Legend.AutoUpdate = 'off';
plot(Handle_Axis1, t_Disturbance_1 * [1, 1], Handle_Axis1.YLim, 'Color', 0.5 * [1, 1, 1], 'LineStyle', '--', 'LineWidth', LineWidth);
plot(Handle_Axis1, t_Disturbance_2 * [1, 1], Handle_Axis1.YLim, 'Color', 0.5 * [1, 1, 1], 'LineStyle', '--', 'LineWidth', LineWidth);

%% Axis for Anti-Windup Variables
Plot_Indeces = [L+3, L+4, L+5, L+6];
Plot_Names = {'$z''_1$', '$z''_2$', '$z''''_1$', '$z''''_2$'};

Handle_Axis2 = axes(Handle_Figure2);
    Handle_Axis2.Box = 'on';
  	Handle_Axis2.BoxStyle = 'full';
 	Handle_Axis2.LineWidth = LineWidth_Thin;
  	Handle_Axis2.FontSize = FontSize;
 	hold(Handle_Axis2, 'on');
 	grid(Handle_Axis2, 'on');
   	Handle_Axis2.XMinorGrid = 'on';
  	Handle_Axis2.YMinorGrid = 'on';
   	Handle_Axis2.Title.String = 'Integral Windup';
  	Handle_Axis2.XTickLabel = [0, t_Disturbance_1, t_Disturbance_2, tf];
  	Handle_Axis2.XLabel.String = 'Time';
    Handle_Axis2.YLabel.String = 'Response';
for i = 1 : length(Plot_Indeces)
    plot(Handle_Axis2, t, X(Plot_Indeces(i), :), 'Color', Colors(i,:), 'LineWidth', LineWidth_Thick);
end
Handle_Legend = legend(Handle_Axis2, Plot_Names);
    Handle_Legend.Interpreter = 'latex';
    Handle_Legend.AutoUpdate = 'off';
plot(Handle_Axis1, t_Disturbance_1 * [1, 1], Handle_Axis1.YLim, 'Color', 0.5 * [1, 1, 1], 'LineStyle', '--', 'LineWidth', LineWidth);
plot(Handle_Axis1, t_Disturbance_2 * [1, 1], Handle_Axis1.YLim, 'Color', 0.5 * [1, 1, 1], 'LineStyle', '--', 'LineWidth', LineWidth);

