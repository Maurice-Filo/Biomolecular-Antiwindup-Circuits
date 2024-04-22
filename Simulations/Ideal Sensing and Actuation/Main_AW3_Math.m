%% Clear Workspace
close all
clear
clc

%% Controller 
% Controller Selection
StoichiometryMatrix_Controller = StoichiometryMatrix_AW3();
PropensityFunction_Controller = @PropensityFunction_AW3;
% Parameters
Parameters_Controller.mu = 10;
Parameters_Controller.eta = 100;
Parameters_Controller.delta_1 = 0;
Parameters_Controller.delta_2 = 0;
M = 2;                  % Number of Controller Species
% Actuation
Parameters_Controller.ActuationParameters.k = 1;
Parameters_Controller.ActuationFunction = @(z_1, Parameters) (Parameters.k * z_1);
% Sensing
Parameters_Controller.SensingParameters.theta = 1;
Parameters_Controller.SensingFunction = @(y, Parameters) (Parameters.theta * y);
% Anti-Windup
Parameters_Controller.WindupParameters.alpha_1 = 1;
Parameters_Controller.WindupParameters.alpha_2 = 1;
Parameters_Controller.WindupParameters.beta_1 = 5;
Parameters_Controller.WindupParameters.beta_2 = 5;
Parameters_Controller.h_1 = @(z_1, Parameters) (1 / (1 + (1/Parameters.alpha_1) * max(z_1 - Parameters.beta_1, 0)));
Parameters_Controller.h_2 = @(z_2, Parameters) (1 / (1 + (1/Parameters.alpha_2) * max(z_2 - Parameters.beta_2, 0)));

%% Reduced Controller
% Controller Selection
StoichiometryMatrix_Controller_Reduced = StoichiometryMatrix_AW3_Reduced();
PropensityFunction_Controller_Reduced = @PropensityFunction_AW3_Reduced; 

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
% Initial Conditions
IC = zeros(L+M,1);      
IC_Reduced = [IC(1:L); IC(L+M) - IC(L+M-1)];    

%% Disturbance 
DisturbanceFactor_1 = 2;
t_Disturbance_1 = 100;
DisturbanceFactor_2 = 0.5;
t_Disturbance_2 = 150;
DisturbedParameter = 'k_0';

%% Ranges
R_Disturbance = [0, Parameters_Plant.gamma_1 * Parameters_Plant.gamma_2 * Parameters_Controller.mu / Parameters_Controller.SensingParameters.theta / Parameters_Plant.k_1]; 
R = [Parameters_Plant.k_1 * Parameters_Plant.k_0 / (Parameters_Plant.gamma_1 * Parameters_Plant.gamma_2), 100];
R_1 = [Parameters_Plant.k_1 * Parameters_Plant.k_0 * DisturbanceFactor_1 / (Parameters_Plant.gamma_1 * Parameters_Plant.gamma_2), 100];
R_2 = [Parameters_Plant.k_1 * Parameters_Plant.k_0 * DisturbanceFactor_1 * DisturbanceFactor_2 / (Parameters_Plant.gamma_1 * Parameters_Plant.gamma_2), 100];

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
                
%% Reduced Closed-Loop Network
% Stoichiometry Matrix
S_Mutual = zeros(L, size(StoichiometryMatrix_Controller_Reduced,2));
S_Mutual(1,Input_Index) = 1;
StoichiometryMatrix_CL_Reduced = [	StoichiometryMatrix_Plant,                                  S_Mutual; ...
                                    zeros(1, size(StoichiometryMatrix_Plant,2)),                StoichiometryMatrix_Controller_Reduced; ...
                                 ];
% Propensity Function
PropensityFunction_CL_Reduced = @(X, Parameters_CL) ...
                                ([ PropensityFunction_Plant(X(1:L), Parameters_CL.Plant); ...
                                PropensityFunction_Controller_Reduced(X(L+1:L+1), X(L), Parameters_CL.Controller); ...
                                ]);
                
%% Simulating the Full & Reduced Models
[t_1, X_1] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, IC, t_Disturbance_1, N_t, Solver);
[~, X_1_Reduced] = DSA(StoichiometryMatrix_CL_Reduced, PropensityFunction_CL_Reduced, Parameters_CL, IC_Reduced, t_Disturbance_1, N_t, Solver);
Parameters_CL.Plant.(DisturbedParameter) = Parameters_CL.Plant.(DisturbedParameter) * DisturbanceFactor_1;
[t_2, X_2] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, X_1(:,end), t_Disturbance_2 - t_Disturbance_1, N_t, Solver);
[~, X_2_Reduced] = DSA(StoichiometryMatrix_CL_Reduced, PropensityFunction_CL_Reduced, Parameters_CL, X_1_Reduced(:,end), t_Disturbance_2 - t_Disturbance_1, N_t, Solver);
Parameters_CL.Plant.(DisturbedParameter) = Parameters_CL.Plant.(DisturbedParameter) * DisturbanceFactor_2;
[t_3, X_3] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, X_2(:,end), tf - t_Disturbance_2, N_t, Solver);
[~, X_3_Reduced] = DSA(StoichiometryMatrix_CL_Reduced, PropensityFunction_CL_Reduced, Parameters_CL, X_2_Reduced(:,end), tf - t_Disturbance_2, N_t, Solver);
t = [t_1, t_Disturbance_1 + t_2(2:end), t_Disturbance_2 + t_3(2:end)];
X = [X_1, X_2(:,2:end), X_3(:,2:end)];
X_Reduced = [X_1_Reduced, X_2_Reduced(:,2:end), X_3_Reduced(:,2:end)];
X_Reduced_Mapped = [X_Reduced(1:L, :); max(X_Reduced(L+1, :), 0); max(-X_Reduced(L+1, :), 0)];

%% Figure Settings
ScalingFactor = 1;
SS = 4;
Figure_Width = 6.5 * SS;
Figure_Height = 4 * SS;
FontSize = ScalingFactor*5 * SS;
FontSize_Small = ScalingFactor*3 * SS;
FontSize_Large = ScalingFactor*6 * SS;
LineWidth = ScalingFactor*0.65 * SS;
LineWidth_Thick = ScalingFactor*0.8 * SS;
LineWidth_Thin = ScalingFactor*0.1 * SS;
MarkerSize = ScalingFactor*2 * SS;
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
    Handle_Figure2.Position = [0, 0, Figure_Width, Figure_Height/3];
    Handle_Figure2.PaperPositionMode = 'auto';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)];

%% Axis for Simulations
Plot_Indeces = [L, L+1, L+2];
Plot_Names = {'$x_L$', '$x_L$ (Reduced)', '$z_1$', '$z_1$ (Reduced)', '$z_2$', '$z_2$ (Reduced)'};
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
  	Handle_Axis1.XTickLabel = {};
    Handle_Axis1.YLabel.String = 'Response';
    Handle_Axis1.YLim = [0, 25];
for i = 1 : length(Plot_Indeces)
    plot(Handle_Axis1, t(1 : floor(length(t)/20) : end), X(Plot_Indeces(i), (1 : floor(length(t)/20) : end)), 'Color', Colors(i,:), 'LineStyle', 'none', 'Marker', 'd', 'MarkerFaceColor', Colors(i,:), 'MarkerSize', MarkerSize);
    plot(Handle_Axis1, t, X_Reduced_Mapped(Plot_Indeces(i), :), 'Color', Colors(i,:), 'LineWidth', LineWidth_Thick, 'LineStyle', '-', 'Marker', 'none');
end
Handle_Legend = legend(Handle_Axis1, Plot_Names);
    Handle_Legend.Interpreter = 'latex';
    Handle_Legend.AutoUpdate = 'off';
Handle_Patch_1 = patch(Handle_Axis1, [0, t_Disturbance_1 * [1, 1], t_Disturbance_2 * [1, 1], tf * [1, 1], ...
                                      t_Disturbance_2 * [1, 1], t_Disturbance_1 * [1, 1], 0, 0], ...
                                      [R(1) * [1, 1], R_1(1) * [1, 1], R_2(1) * [1, 1], ...
                                       R_2(2) * [1, 1], R_1(2) * [1, 1], R(2) * [1, 1], R(1)], ...
                                      Colors(1,:), 'LineStyle', '-');
    Handle_Patch_1.FaceAlpha = 0.2;
plot(Handle_Axis1, Handle_Axis1.XLim, Parameters_Controller.mu / Parameters_Controller.SensingParameters.theta * [1, 1], 'LineWidth', LineWidth, 'Color', 'k', 'LineStyle', '--');

%% Axis for Input/Disturbance
Handle_Axis2 = axes(Handle_Figure2);
    Handle_Axis2.Box = 'on';
  	Handle_Axis2.BoxStyle = 'full';
 	Handle_Axis2.LineWidth = LineWidth_Thin;
  	Handle_Axis2.FontSize = FontSize;
 	hold(Handle_Axis2, 'on');
 	grid(Handle_Axis2, 'on');
   	Handle_Axis2.XMinorGrid = 'on';
  	Handle_Axis2.YMinorGrid = 'on';
  	Handle_Axis2.XLabel.String = 'Time';
    Handle_Axis2.YLabel.String = '$\Delta(t)$';
    Handle_Axis2.YLabel.Interpreter = 'latex';
Handle_Axis2.YLim = [0*min([Parameters_Plant.(DisturbedParameter), Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1, Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1 * DisturbanceFactor_2]), ...
                     1.5*max([Parameters_Plant.(DisturbedParameter), Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1, Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1 * DisturbanceFactor_2])];
plot(Handle_Axis2, [0, t_Disturbance_1 * [1, 1], t_Disturbance_2 * [1, 1], tf], ...
                   [Parameters_Plant.(DisturbedParameter) * [1, 1], ...
                    Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1 * [1, 1], ...
                    Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1 * DisturbanceFactor_2 * [1, 1]], ...
                    'LineWidth', LineWidth, 'Color', 0.5*[1, 1, 1]);
Handle_Patch_2 = patch(Handle_Axis2, [Handle_Axis2.XLim, flip(Handle_Axis2.XLim)], [R_Disturbance(1), R_Disturbance(1), R_Disturbance(2), R_Disturbance(2)], 0.5 * [1, 1, 1], 'LineStyle', '-');
    Handle_Patch_2.FaceAlpha = 0.2;
    
%% Save Figures
Save_Flag = 0;
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
    print(Handle_Figure1, 'Simulation_NEC', '-dpdf', '-painters');
    Handle_Figure1.Color = [1, 1, 1];
    
    Handle_Figure2.Color = 'none';
    set(Handle_Figure2, 'InvertHardCopy', 'off');
    print(Handle_Figure2, 'Simulation_NEC_Disturbance', '-dpdf', '-painters');
    Handle_Figure2.Color = [1, 1, 1];
end
