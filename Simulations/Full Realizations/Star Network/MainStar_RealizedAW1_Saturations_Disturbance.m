%% Clear Workspace
close all
clear
clc

Factor = 1e3;

%% Controller 
% Controller Selection
StoichiometryMatrix_Controller = StoichiometryMatrix_PID3_AW1_Realization();
PropensityFunction_Controller = @PropensityFunction_PID3_AW1_Realization;
% Parameters
Parameters_Controller.mu = 6;
Parameters_Controller.eta = 1e4;
Parameters_Controller.delta_1 = 0;
Parameters_Controller.delta_2 = 0;
Parameters_Controller.kappa = 1e-4;
Parameters_Controller.kappa_0 = 10;
Parameters_Controller.alpha_0 = 1;
Parameters_Controller.gamma_0 = 0.2;
Parameters_Controller.delta_0 = 1.6;
Parameters_Controller.delta = 0.01;
M = 7;                  % Number of Controller Species
% Actuation
Parameters_Controller.ActuationParameters.k = 0.015 * Factor;
Parameters_Controller.ActuationParameters.kappa_a = Factor;
Parameters_Controller.ActuationFunction = @(z_1, Parameters) (Parameters.k * z_1 / (z_1 + Parameters.kappa_a));
% Sensing
Parameters_Controller.SensingParameters.theta = 0.3 * Factor;
Parameters_Controller.SensingParameters.kappa_s = Factor;
Parameters_Controller.SensingFunction = @(y, Parameters) (Parameters.theta * y / (y + Parameters.kappa_s));
% Anti-Windup
Parameters_Controller.eta_v = 1e4;
Parameters_Controller.eta_w = 1e4;
Parameters_Controller.v_0 = 3500;
Parameters_Controller.w_0 = 3500;
Parameters_Controller.alpha_v = 1;
Parameters_Controller.alpha_w = 1;
Parameters_Controller.delta_v1 = 1;
Parameters_Controller.delta_v2 = 1;
Parameters_Controller.delta_w1 = 1;
Parameters_Controller.delta_w2 = 1;
Parameters_Controller.h_1 = @(v_1, v_2) Parameters_Controller.alpha_v*v_2;
Parameters_Controller.h_2 = @(w_1, w_2) Parameters_Controller.alpha_w*w_2;
Parameters_Controller.g_1 = @(z_2) z_2;
Parameters_Controller.g_2 = @(z_1) z_1;
% Without Anti-Windup
Parameters_Controller_NoAW = Parameters_Controller;
Parameters_Controller_NoAW.alpha_v = 0;
Parameters_Controller_NoAW.alpha_w = 0;
Parameters_Controller_NoAW.h_1 = @(v_1, v_2) Parameters_Controller_NoAW.alpha_v*v_2;
Parameters_Controller_NoAW.h_2 = @(w_1, w_2) Parameters_Controller_NoAW.alpha_w*w_2;

%% Plant
% Plant Selection
StoichiometryMatrix_Plant = StoichiometryMatrix_Star();
PropensityFunction_Plant = @PropensityFunction_Star;
% Parameters
Parameters_Plant.b_0 = .2;
Parameters_Plant.b_1 = .2;
Parameters_Plant.b_2 = .2;
Parameters_Plant.b_3 = .2;
Parameters_Plant.b_4 = .2;
Parameters_Plant.b_5 = .2;
Parameters_Plant.k_1 = 0.1 * Factor;
Parameters_Plant.k_2 = 0.1 * Factor;
Parameters_Plant.k_3 = 0.1 * Factor;
Parameters_Plant.k_4 = 0.1 * Factor;
Parameters_Plant.k_5 = 0.1 * Factor;
Parameters_Plant.kappa_1 = Factor;
Parameters_Plant.kappa_2 = Factor;
Parameters_Plant.kappa_3 = Factor;
Parameters_Plant.kappa_4 = Factor;
Parameters_Plant.kappa_5 = Factor;
Parameters_Plant.gamma_F = 0.3;
Parameters_Plant.kappa_F = 1e-2;
Parameters_Plant.gamma_1 = 0.1;
Parameters_Plant.gamma_2 = 0.1;
Parameters_Plant.gamma_3 = 0.1;
Parameters_Plant.gamma_4 = 0.1;
Parameters_Plant.gamma_5 = 0.1;
Parameters_Plant.gamma_6 = 0.1;
L = 6;                  % Number of Plant Species
Input_Index = 1;
Output_Index = 6;

%% Simulation Settings
tf = 20000;             % Final Time
N_t = 10000;             % Time Samples
Solver = 'ODE15s';      % ODE Solver
IC = zeros(L+M,1);      % Initial Condition  

%% Disturbance 
DisturbanceFactor_1 = 1.3;
t_Disturbance_1 = 5000;

DisturbanceFactor_2 = 0.1;
t_Disturbance_2 = 10000;

DisturbanceFactor_3 = 10;
t_Disturbance_3 = 11000;

DisturbanceFactor_4 = 1;
t_Disturbance_4 = 19000;

DisturbedParameter0 = 'gamma_6';
DisturbedParameter1 = 'k_1';
DisturbedParameter2 = 'k_2';
DisturbedParameter3 = 'k_3';
DisturbedParameter4 = 'k_4';
DisturbedParameter5 = 'k_5';

%% Ranges
r_in = Parameters_Controller.mu / Parameters_Controller.SensingParameters.theta;
r_out = Parameters_Controller.SensingParameters.kappa_s * r_in / (1 - r_in);
               
%% Closed-Loop Network with Anti-Windup
% Parameters                    
Parameters_CL.Controller = Parameters_Controller;
Parameters_CL.Plant = Parameters_Plant;
% Stoichiometry Matrix
S_Mutual = zeros(L, size(StoichiometryMatrix_Controller,2));
S_Mutual(Input_Index, 1) = 1;
S_Mutual(Input_Index,end) = -1;
StoichiometryMatrix_CL = [	StoichiometryMatrix_Plant,                          S_Mutual; ...
                            zeros(M, size(StoichiometryMatrix_Plant,2)),        StoichiometryMatrix_Controller; ...
                         ];
% Propensity Function
PropensityFunction_CL = @(X, Parameters_CL) ...
                    ([ PropensityFunction_Plant(X(1:L), Parameters_CL.Plant); ...
                       PropensityFunction_Controller(X(L+1:L+M), X(L), X(1), Parameters_CL.Controller); ...
                    ]);  
                
%% Closed-Loop Network without Anti-Windup
% Parameters                    
Parameters_CL_NoAW.Controller = Parameters_Controller_NoAW;
Parameters_CL_NoAW.Plant = Parameters_Plant;

%% Simulating the Full Model with Anti-Windup
[t_1, X_1] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, IC, t_Disturbance_1, N_t, Solver);

Parameters_CL.Plant.(DisturbedParameter0) = Parameters_CL.Plant.(DisturbedParameter0) * DisturbanceFactor_1;
[t_2, X_2] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, X_1(:,end), t_Disturbance_2 - t_Disturbance_1, N_t, Solver);

Parameters_CL.Plant.(DisturbedParameter1) = Parameters_CL.Plant.(DisturbedParameter1) * DisturbanceFactor_2;
Parameters_CL.Plant.(DisturbedParameter2) = Parameters_CL.Plant.(DisturbedParameter2) * DisturbanceFactor_2;
Parameters_CL.Plant.(DisturbedParameter3) = Parameters_CL.Plant.(DisturbedParameter3) * DisturbanceFactor_2;
Parameters_CL.Plant.(DisturbedParameter4) = Parameters_CL.Plant.(DisturbedParameter4) * DisturbanceFactor_2;
Parameters_CL.Plant.(DisturbedParameter5) = Parameters_CL.Plant.(DisturbedParameter5) * DisturbanceFactor_2;
[t_3, X_3] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, X_2(:,end), t_Disturbance_3 - t_Disturbance_2, N_t, Solver);

Parameters_CL.Plant.(DisturbedParameter1) = Parameters_CL.Plant.(DisturbedParameter1) * DisturbanceFactor_3;
Parameters_CL.Plant.(DisturbedParameter2) = Parameters_CL.Plant.(DisturbedParameter2) * DisturbanceFactor_3;
Parameters_CL.Plant.(DisturbedParameter3) = Parameters_CL.Plant.(DisturbedParameter3) * DisturbanceFactor_3;
Parameters_CL.Plant.(DisturbedParameter4) = Parameters_CL.Plant.(DisturbedParameter4) * DisturbanceFactor_3;
Parameters_CL.Plant.(DisturbedParameter5) = Parameters_CL.Plant.(DisturbedParameter5) * DisturbanceFactor_3;
[t_4, X_4] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, X_3(:,end), t_Disturbance_4 - t_Disturbance_3, N_t, Solver);

Parameters_CL.Plant.(DisturbedParameter1) = Parameters_CL.Plant.(DisturbedParameter1) * DisturbanceFactor_4;
[t_5, X_5] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, X_4(:,end), tf - t_Disturbance_4, N_t, Solver);

t_AW = [t_1, t_Disturbance_1 + t_2(2:end), t_Disturbance_2 + t_3(2:end), t_Disturbance_3 + t_4(2:end), t_Disturbance_4 + t_5(2:end)];
X_AW = [X_1, X_2(:,2:end), X_3(:,2:end), X_4(:,2:end), X_5(:,2:end)];

%% Simulating the Full Model without Anti-Windup
[t_1, X_1] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL_NoAW, IC, t_Disturbance_1, N_t, Solver);

Parameters_CL_NoAW.Plant.(DisturbedParameter0) = Parameters_CL_NoAW.Plant.(DisturbedParameter0) * DisturbanceFactor_1;
[t_2, X_2] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL_NoAW, X_1(:,end), t_Disturbance_2 - t_Disturbance_1, N_t, Solver);

Parameters_CL_NoAW.Plant.(DisturbedParameter1) = Parameters_CL_NoAW.Plant.(DisturbedParameter1) * DisturbanceFactor_2;
Parameters_CL_NoAW.Plant.(DisturbedParameter2) = Parameters_CL_NoAW.Plant.(DisturbedParameter2) * DisturbanceFactor_2;
Parameters_CL_NoAW.Plant.(DisturbedParameter3) = Parameters_CL_NoAW.Plant.(DisturbedParameter3) * DisturbanceFactor_2;
Parameters_CL_NoAW.Plant.(DisturbedParameter4) = Parameters_CL_NoAW.Plant.(DisturbedParameter4) * DisturbanceFactor_2;
Parameters_CL_NoAW.Plant.(DisturbedParameter5) = Parameters_CL_NoAW.Plant.(DisturbedParameter5) * DisturbanceFactor_2;
[t_3, X_3] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL_NoAW, X_2(:,end), t_Disturbance_3 - t_Disturbance_2, N_t, Solver);

Parameters_CL_NoAW.Plant.(DisturbedParameter1) = Parameters_CL_NoAW.Plant.(DisturbedParameter1) * DisturbanceFactor_3;
Parameters_CL_NoAW.Plant.(DisturbedParameter2) = Parameters_CL_NoAW.Plant.(DisturbedParameter2) * DisturbanceFactor_3;
Parameters_CL_NoAW.Plant.(DisturbedParameter3) = Parameters_CL_NoAW.Plant.(DisturbedParameter3) * DisturbanceFactor_3;
Parameters_CL_NoAW.Plant.(DisturbedParameter4) = Parameters_CL_NoAW.Plant.(DisturbedParameter4) * DisturbanceFactor_3;
Parameters_CL_NoAW.Plant.(DisturbedParameter5) = Parameters_CL_NoAW.Plant.(DisturbedParameter5) * DisturbanceFactor_3;
[t_4, X_4] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL_NoAW, X_3(:,end), t_Disturbance_4 - t_Disturbance_3, N_t, Solver);

Parameters_CL_NoAW.Plant.(DisturbedParameter1) = Parameters_CL_NoAW.Plant.(DisturbedParameter1) * DisturbanceFactor_4;
[t_5, X_5] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL_NoAW, X_4(:,end), tf - t_Disturbance_4, N_t, Solver);

t = [t_1, t_Disturbance_1 + t_2(2:end), t_Disturbance_2 + t_3(2:end), t_Disturbance_3 + t_4(2:end), t_Disturbance_4 + t_5(2:end)];
X = [X_1, X_2(:,2:end), X_3(:,2:end), X_4(:,2:end), X_5(:,2:end)];

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
Colors = [ ...
          0, 152, 58; ... % Green
          39, 52, 139; ... % Blue
          227, 6, 19; ... % Red
          230, 0, 126; ... % Magenta
          130, 54, 140; ... % Purple
          ]/255;
    
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
    Handle_Figure2.Position = [0, 0, Figure_Width, Figure_Height/2.5];
    Handle_Figure2.PaperPositionMode = 'auto';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)];
    
%% Set Figure 3
% Handle_Figure3 = figure();
%     Handle_Figure3.Color = [1 1 1];
%     Handle_Figure3.PaperUnits = 'centimeters';
%     Handle_Figure3.Units = 'centimeters';
%     Handle_Figure3.Position = [0, 0, Figure_Width, Figure_Height * 1.5];
%     Handle_Figure3.PaperPositionMode = 'auto';
%     Handle_Figure3.PaperSize = [Handle_Figure3.PaperPosition(3), Handle_Figure3.PaperPosition(4)];

%% Axis for Simulations
Plot_Indeces = [L, L+1, L+2];
Plot_Names = {'x_6 without Anti-Windup', 'x_6 with Anti-Windup', 'z_1 without Anti-Windup', 'z_1 with Anti-Windup'};

Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.09, 0.05, 0.815, 0.92];
    Handle_Axis1.Box = 'on';
  	Handle_Axis1.BoxStyle = 'full';
 	Handle_Axis1.LineWidth = LineWidth_Thin;
  	Handle_Axis1.FontSize = FontSize;
 	hold(Handle_Axis1, 'on');
 	grid(Handle_Axis1, 'on');
   	Handle_Axis1.XMinorGrid = 'on';
  	Handle_Axis1.YMinorGrid = 'on';
  	Handle_Axis1.XTickLabel = {};
    Handle_Axis1.YLabel.String = 'Response';
    Handle_Axis1.YLim = [0, 40];
yyaxis(Handle_Axis1, 'left');
Handle_Axis1.YLim = [0, 30];
plot(Handle_Axis1, t, X(Plot_Indeces(1), :), 'LineStyle', '-', 'Color', [Colors(1,:), 0.3], 'LineWidth', LineWidth_Thick*2, 'Marker', 'none');
plot(Handle_Axis1, t_AW, X_AW(Plot_Indeces(1), :), 'Color', Colors(1,:), 'LineWidth', LineWidth_Thick/1.5, 'LineStyle', '-', 'Marker', 'none');
yyaxis(Handle_Axis1, 'right'); 
Handle_Axis1.YAxis(1).Color = 'k';
Handle_Axis1.YAxis(2).Color = Colors(2,:);
plot(Handle_Axis1, t, X(Plot_Indeces(2), :), 'LineStyle', '-', 'Color', [Colors(2,:), 0.3], 'LineWidth', LineWidth_Thick*2, 'LineStyle', '-', 'Marker', 'none');
plot(Handle_Axis1, t_AW, X_AW(Plot_Indeces(2), :), 'Color', Colors(2,:), 'LineWidth', LineWidth_Thick/1.5, 'LineStyle', '-', 'Marker', 'none');
yyaxis(Handle_Axis1, 'left');
% plot(Handle_Axis1, t_AW, X_AW(Plot_Indeces(3), :), 'Color', Colors(3,:), 'LineWidth', LineWidth_Thick, 'LineStyle', '-', 'Marker', 'none');

Handle_Legend = legend(Handle_Axis1, Plot_Names);
%     Handle_Legend.Interpreter = 'latex';
    Handle_Legend.AutoUpdate = 'off';
    Handle_Legend.Position = [0.63, 0.1, 0.2, 0.2];
    
plot(Handle_Axis1, [0, tf], r_out * [1, 1], 'LineWidth', LineWidth', 'Color', 'k', 'LineStyle', '--');

% Inset
Handle_Axis11 = axes(Handle_Figure1);
    Handle_Axis11.Position = [0.12, 0.77, 0.37, 0.14];
    Handle_Axis11.Box = 'on';
  	Handle_Axis11.BoxStyle = 'full';
 	Handle_Axis11.LineWidth = LineWidth_Thin;
  	Handle_Axis11.FontSize = FontSize/1.2;
 	hold(Handle_Axis11, 'on');
 	grid(Handle_Axis11, 'on');
   	Handle_Axis11.XMinorGrid = 'on';
  	Handle_Axis11.YMinorGrid = 'on';
    Handle_Axis11.XTick = [0, 0.5, 1, 1.5, 2] * 1e4;
  	Handle_Axis11.XTickLabel = {0, 0.5, 1, 1.5, '2x10^4'};
    Handle_Axis11.Title.String = 'Anti-Windup Species W_2';
    Handle_Axis11.Title.Position(1:2) = [1e4, 6.7];
plot(Handle_Axis11, t_AW, X_AW(L+6,:), 'Color', Colors(5,:), 'LineWidth', LineWidth);
    
%% Axis for Input/Disturbance
Handle_Axis2 = axes(Handle_Figure2);
    Handle_Axis2.Position = [0.09, 0.35, 0.815, 0.6];
    Handle_Axis2.Box = 'on';
  	Handle_Axis2.BoxStyle = 'full';
 	Handle_Axis2.LineWidth = LineWidth_Thin;
  	Handle_Axis2.FontSize = FontSize;
 	hold(Handle_Axis2, 'on');
 	grid(Handle_Axis2, 'on');
   	Handle_Axis2.XMinorGrid = 'off';
  	Handle_Axis2.YMinorGrid = 'off';
  	Handle_Axis2.XLabel.String = 'Time';
%     Handle_Axis2.YTick = [2, 5, 8, 11];
% Handle_Axis2.YLim = [0.9*min([Parameters_Plant.(DisturbedParameter1), Parameters_Plant.(DisturbedParameter1) * DisturbanceFactor_1, Parameters_Plant.(DisturbedParameter1) * DisturbanceFactor_1 * DisturbanceFactor_2, Parameters_Plant.(DisturbedParameter1) * DisturbanceFactor_1 * DisturbanceFactor_2 * DisturbanceFactor_3, Parameters_Plant.(DisturbedParameter1) * DisturbanceFactor_1 * DisturbanceFactor_2 * DisturbanceFactor_3]), ...
%                      1.1*max([Parameters_Plant.(DisturbedParameter1), Parameters_Plant.(DisturbedParameter1) * DisturbanceFactor_1, Parameters_Plant.(DisturbedParameter1) * DisturbanceFactor_1 * DisturbanceFactor_2, Parameters_Plant.(DisturbedParameter1) * DisturbanceFactor_1 * DisturbanceFactor_2 * DisturbanceFactor_3, Parameters_Plant.(DisturbedParameter1) * DisturbanceFactor_1 * DisturbanceFactor_2 * DisturbanceFactor_3])];
yyaxis(Handle_Axis2, 'left');
Handle_Axis2.YAxis(1).Color = [0, 0, 0];
Handle_Axis2.YTick = [0.05, 0.1, 0.15];
Handle_Axis2.YLim = [0.05, 0.15];
Handle_Axis2.YLabel.String = '$\gamma_6(t)$';
Handle_Axis2.YLabel.Position(1) = -0.13e4;
Handle_Axis2.YLabel.Interpreter = 'latex';
plot(Handle_Axis2, [0, t_Disturbance_1 * [1, 1], tf], [Parameters_Plant.(DisturbedParameter0) * [1, 1], Parameters_Plant.(DisturbedParameter0) * [1, 1] * DisturbanceFactor_1], 'LineWidth', LineWidth, 'Color', 'k');
yyaxis(Handle_Axis2, 'right');
Handle_Axis2.YAxis(2).Color = Colors(4,:);
Handle_Axis2.YTick = [0, 50, 100];
Handle_Axis2.YLim = [0, 110];
Handle_Axis2.YLabel.String = '$k_p(t)$';
Handle_Axis2.YLabel.Interpreter = 'latex';
plot(Handle_Axis2, [0, t_Disturbance_2 * [1, 1], t_Disturbance_3 * [1, 1], tf], ...
                   [Parameters_Plant.(DisturbedParameter1) * [1, 1], ...
                    Parameters_Plant.(DisturbedParameter1) * DisturbanceFactor_2 * [1, 1], ...
                    Parameters_Plant.(DisturbedParameter1) * DisturbanceFactor_2 * DisturbanceFactor_3 * [1, 1]], ...
                    'LineWidth', LineWidth, 'Color', Colors(4,:));
    
%% Save Figures
Save_Flag = 1;
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
    print(Handle_Figure1, 'PID3_RealizedAW1_Saturation_Rejection', '-dpdf', '-painters');
    Handle_Figure1.Color = [1, 1, 1];
    
    Handle_Figure2.Color = 'none';
    set(Handle_Figure2, 'InvertHardCopy', 'off');
    print(Handle_Figure2, 'PID3_RealizedAW1_Saturation_Rejection_Disturbance1', '-dpdf', '-painters');
    Handle_Figure2.Color = [1, 1, 1];
end
