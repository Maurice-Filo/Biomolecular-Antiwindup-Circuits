%% Clear Workspace
close all
clear
clc

%% Controller 
% Controller Selection
StoichiometryMatrix_Controller = StoichiometryMatrix_AW2_Realization();
PropensityFunction_Controller = @PropensityFunction_AW2_Realization;
% Parameters
Parameters_Controller.mu = 10;
Parameters_Controller.eta = 100;
Parameters_Controller.delta_1 = 0;
Parameters_Controller.delta_2 = 0;
M = 6;                  % Number of Controller Species
% Actuation
Parameters_Controller.ActuationParameters.k = 8;
Parameters_Controller.ActuationParameters.kappa_a = 5;
Parameters_Controller.ActuationFunction = @(z_1, Parameters) (Parameters.k * z_1 / (z_1 + Parameters.kappa_a));
% Sensing
Parameters_Controller.SensingParameters.theta = 15;
Parameters_Controller.SensingParameters.kappa_s = 5;
Parameters_Controller.SensingFunction = @(y, Parameters) (Parameters.theta * y / (y + Parameters.kappa_s));
% Anti-Windup
Parameters_Controller.eta_v = 100;
Parameters_Controller.eta_w = 100;
Parameters_Controller.v_0 = 20;
Parameters_Controller.w_0 = 10;
Parameters_Controller.alpha_v = 1;
Parameters_Controller.alpha_w = 1;
Parameters_Controller.delta_v1 = 1;
Parameters_Controller.delta_v2 = 1;
Parameters_Controller.delta_w1 = 1;
Parameters_Controller.delta_w2 = 1;
Parameters_Controller.h_1 = @(v_1, v_2) Parameters_Controller.alpha_v*v_2;
Parameters_Controller.h_2 = @(w_1, w_2) Parameters_Controller.alpha_w*w_2;
Parameters_Controller.g_1 = @(z_1) z_1;
Parameters_Controller.g_2 = @(z_2) z_2;

%% Plant
% Plant Selection
StoichiometryMatrix_Plant = StoichiometryMatrix_GeneExp();
PropensityFunction_Plant = @PropensityFunction_GeneExp;
% Parameters
Parameters_Plant.k_0 = 5;
Parameters_Plant.k_1 = 1;
Parameters_Plant.gamma_1 = 1;
Parameters_Plant.gamma_2 = 1;
L = 2;                  % Number of Plant Species
Input_Index = 1;
Output_Index = 2;

%% Simulation Settings
tf = 500;               % Final Time
N_t = 1000;             % Time Samples
Solver = 'ODE23s';      % ODE Solver
IC = zeros(L+M,1);      % Initial Condition   

%% Disturbance 
DisturbanceFactor_1 = 2.25;
t_Disturbance_1 = 50;
DisturbanceFactor_2 = 1/2.25;
t_Disturbance_2 = 175;
DisturbanceFactor_3 = 0.1;
t_Disturbance_3 = 300;
DisturbanceFactor_4 = 10;
t_Disturbance_4 = 350;
DisturbedParameter = 'k_0';

%% Ranges
r_in = Parameters_Controller.mu / Parameters_Controller.SensingParameters.theta;
r_out = Parameters_Controller.SensingParameters.kappa_s * r_in / (1 - r_in);
R_Disturbance = [max(Parameters_Plant.gamma_1*Parameters_Plant.gamma_2*r_out/Parameters_Plant.k_1 - Parameters_Controller.ActuationParameters.k, 0), ...
                 Parameters_Plant.gamma_1*Parameters_Plant.gamma_2*r_out/Parameters_Plant.k_1];
R = [Parameters_Plant.k_1 * Parameters_Plant.k_0 / (Parameters_Plant.gamma_1 * Parameters_Plant.gamma_2), Parameters_Plant.k_1 * (Parameters_Plant.k_0 + Parameters_Controller.ActuationParameters.k) / (Parameters_Plant.gamma_1 * Parameters_Plant.gamma_2)];
R_1 = [Parameters_Plant.k_1 * Parameters_Plant.k_0 * DisturbanceFactor_1 / (Parameters_Plant.gamma_1 * Parameters_Plant.gamma_2), Parameters_Plant.k_1 * (Parameters_Plant.k_0 * DisturbanceFactor_1 + Parameters_Controller.ActuationParameters.k) / (Parameters_Plant.gamma_1 * Parameters_Plant.gamma_2)];
R_2 = [Parameters_Plant.k_1 * Parameters_Plant.k_0 * DisturbanceFactor_1 * DisturbanceFactor_2 / (Parameters_Plant.gamma_1 * Parameters_Plant.gamma_2), Parameters_Plant.k_1 * (Parameters_Plant.k_0 * DisturbanceFactor_1 * DisturbanceFactor_2 + Parameters_Controller.ActuationParameters.k) / (Parameters_Plant.gamma_1 * Parameters_Plant.gamma_2)];
R_3 = [Parameters_Plant.k_1 * Parameters_Plant.k_0 * DisturbanceFactor_1 * DisturbanceFactor_2 * DisturbanceFactor_3 / (Parameters_Plant.gamma_1 * Parameters_Plant.gamma_2), Parameters_Plant.k_1 * (Parameters_Plant.k_0 * DisturbanceFactor_1 * DisturbanceFactor_2 * DisturbanceFactor_3 + Parameters_Controller.ActuationParameters.k) / (Parameters_Plant.gamma_1 * Parameters_Plant.gamma_2)];
R_4 = [Parameters_Plant.k_1 * Parameters_Plant.k_0 * DisturbanceFactor_1 * DisturbanceFactor_2 * DisturbanceFactor_3 * DisturbanceFactor_4 / (Parameters_Plant.gamma_1 * Parameters_Plant.gamma_2), Parameters_Plant.k_1 * (Parameters_Plant.k_0 * DisturbanceFactor_1 * DisturbanceFactor_2 * DisturbanceFactor_3 * DisturbanceFactor_4 + Parameters_Controller.ActuationParameters.k) / (Parameters_Plant.gamma_1 * Parameters_Plant.gamma_2)];
               
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
[t_3, X_3] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, X_2(:,end), t_Disturbance_3 - t_Disturbance_2, N_t, Solver);

Parameters_CL.Plant.(DisturbedParameter) = Parameters_CL.Plant.(DisturbedParameter) * DisturbanceFactor_3;
[t_4, X_4] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, X_3(:,end), t_Disturbance_4 - t_Disturbance_3, N_t, Solver);

Parameters_CL.Plant.(DisturbedParameter) = Parameters_CL.Plant.(DisturbedParameter) * DisturbanceFactor_4;
[t_5, X_5] = DSA(StoichiometryMatrix_CL, PropensityFunction_CL, Parameters_CL, X_4(:,end), tf - t_Disturbance_4, N_t, Solver);

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
    Handle_Figure2.Position = [0, 0, Figure_Width, Figure_Height/3];
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
Plot_Names = {'$x_L$', '$z_1$', '$z_2$'};

Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Position = [0.08, 0.05, 0.89, 0.92];
    Handle_Axis1.Box = 'on';
  	Handle_Axis1.BoxStyle = 'full';
 	Handle_Axis1.LineWidth = LineWidth_Thin;
  	Handle_Axis1.FontSize = FontSize;
 	hold(Handle_Axis1, 'on');
 	grid(Handle_Axis1, 'on');
   	Handle_Axis1.XMinorGrid = 'on';
  	Handle_Axis1.YMinorGrid = 'on';
%    	Handle_Axis1.Title.String = 'Disturbance Rejection: Realized Anti-Windup Topology I';
  	Handle_Axis1.XTickLabel = {};
    Handle_Axis1.YLabel.String = 'Response';
    Handle_Axis1.YLim = [0, 25];
for i = 1 : length(Plot_Indeces)
    plot(Handle_Axis1, t, X(Plot_Indeces(i), :), 'Color', Colors(i,:), 'LineWidth', LineWidth_Thick, 'LineStyle', '-', 'Marker', 'none');
end
Handle_Legend = legend(Handle_Axis1, Plot_Names);
    Handle_Legend.Interpreter = 'latex';
    Handle_Legend.AutoUpdate = 'off';
    Handle_Legend.Location = 'best';
    
Handle_Patch_1 = patch(Handle_Axis1, [0, t_Disturbance_1 * [1, 1], t_Disturbance_2 * [1, 1], t_Disturbance_3 * [1, 1], t_Disturbance_4 * [1, 1], tf, ...
                                      tf, t_Disturbance_4 * [1, 1],  t_Disturbance_3 * [1, 1],  t_Disturbance_2 * [1, 1],  t_Disturbance_1 * [1, 1], 0], ...
                                     [R(1) * [1, 1], R_1(1) * [1, 1], R_2(1) * [1, 1], R_3(1) * [1, 1], R_4(1) * [1, 1], ...
                                      R_4(2) * [1, 1], R_3(2) * [1, 1], R_2(2) * [1, 1], R_1(2) * [1, 1], R(2) * [1, 1]], Colors(1,:) , 'LineStyle', '-');
    Handle_Patch_1.FaceAlpha = 0.2;
    
plot(Handle_Axis1, [0, tf], r_out * [1, 1], 'LineWidth', LineWidth', 'Color', 'k', 'LineStyle', '--');

% Inset
Handle_Axis11 = axes(Handle_Figure1);
    Handle_Axis11.Position = [0.12, 0.77, 0.48, 0.14];
    Handle_Axis11.Box = 'on';
  	Handle_Axis11.BoxStyle = 'full';
 	Handle_Axis11.LineWidth = LineWidth_Thin;
  	Handle_Axis11.FontSize = FontSize/1.2;
 	hold(Handle_Axis11, 'on');
 	grid(Handle_Axis11, 'on');
   	Handle_Axis11.XMinorGrid = 'on';
  	Handle_Axis11.YMinorGrid = 'on';
    Handle_Axis11.Title.String = 'Anti-Windup Species: V_2 and W_2';
    Handle_Axis11.Title.Position(1:2) = [250, .1];
plot(Handle_Axis11, t, X(L+4,:), 'Color', Colors(5,:), 'LineWidth', LineWidth);
plot(Handle_Axis11, t, X(L+6,:), 'Color', [Colors(5,:), 0.5], 'LineWidth', LineWidth);
Handle_Legend_Inset = legend(Handle_Axis11, {'V_2', 'W_2'});
    
%% Axis for Input/Disturbance
Handle_Axis2 = axes(Handle_Figure2);
    Handle_Axis2.Position = [0.08, 0.35, 0.89, 0.6];
    Handle_Axis2.Box = 'on';
  	Handle_Axis2.BoxStyle = 'full';
 	Handle_Axis2.LineWidth = LineWidth_Thin;
  	Handle_Axis2.FontSize = FontSize;
 	hold(Handle_Axis2, 'on');
 	grid(Handle_Axis2, 'on');
   	Handle_Axis2.XMinorGrid = 'off';
  	Handle_Axis2.YMinorGrid = 'off';
  	Handle_Axis2.XLabel.String = 'Time';
    Handle_Axis2.YLabel.String = '$\Delta(t)$';
    Handle_Axis2.YLabel.Interpreter = 'latex';
    Handle_Axis2.YTick = [2, 5, 8, 11];
Handle_Axis2.YLim = [0.9*min([Parameters_Plant.(DisturbedParameter), Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1, Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1 * DisturbanceFactor_2, Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1 * DisturbanceFactor_2 * DisturbanceFactor_3, Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1 * DisturbanceFactor_2 * DisturbanceFactor_3]), ...
                     1.1*max([Parameters_Plant.(DisturbedParameter), Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1, Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1 * DisturbanceFactor_2, Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1 * DisturbanceFactor_2 * DisturbanceFactor_3, Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1 * DisturbanceFactor_2 * DisturbanceFactor_3])];
plot(Handle_Axis2, [0, t_Disturbance_1 * [1, 1], t_Disturbance_2 * [1, 1], t_Disturbance_3 * [1, 1], t_Disturbance_4 * [1, 1] tf], ...
                   [Parameters_Plant.(DisturbedParameter) * [1, 1], ...
                    Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1 * [1, 1], ...
                    Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1 * DisturbanceFactor_2 * [1, 1], ...
                    Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1 * DisturbanceFactor_2 * DisturbanceFactor_3 * [1, 1], ...
                    Parameters_Plant.(DisturbedParameter) * DisturbanceFactor_1 * DisturbanceFactor_2 * DisturbanceFactor_3 * DisturbanceFactor_4 * [1, 1]], ...
                    'LineWidth', LineWidth, 'Color', Colors(4,:));
                
Handle_Patch_2 = patch(Handle_Axis2, [0, tf, tf, 0], [R_Disturbance(1)*[1, 1], R_Disturbance(2)*[1, 1]], 0.5*[1, 1, 1], 'LineStyle', '-');
    Handle_Patch_2.FaceAlpha = 0.2;
    
%% Axis for Simulations of Anti-Windup Variables
% Plot_Indeces_AW = [L+3, L+4, L+5, L+6];
% Plot_Names_AW = {'$v_1$', '$v_2$', '$w_1$', '$w_2$'};
% 
% figure(Handle_Figure3);
% Handle_Axis3 = gobjects(4,1);
% for j = 1 : 4
% Handle_Axis3(j) = subplot(4, 1, j);
%     Handle_Axis3(j).Box = 'on';
%   	Handle_Axis3(j).BoxStyle = 'full';
%  	Handle_Axis3(j).LineWidth = LineWidth_Thin;
%   	Handle_Axis3(j).FontSize = FontSize;
%  	hold(Handle_Axis3(j), 'on');
%  	grid(Handle_Axis3(j), 'on');
%    	Handle_Axis3(j).XMinorGrid = 'on';
%   	Handle_Axis3(j).YMinorGrid = 'on';
%     if j~= 4
%         Handle_Axis3(j).XTickLabel = {};
%     end
%     Handle_Axis3(j).YLabel.String = Plot_Names_AW{j};
%     Handle_Axis3(j).YLabel.Interpreter = 'latex';
% plot(Handle_Axis3(j), t, X(Plot_Indeces_AW(j), :), 'Color', Colors(5,:), 'LineWidth', LineWidth_Thick);
% end
% Handle_Axis3(1).Title.String = 'Topology II: Anti-Windup Species';
% Handle_Axis3(4).XLabel.String = 'Time';
    
%% Save Figures
Save_Flag = 1;
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
    print(Handle_Figure1, 'RealizedAW2_Saturation_Rejection', '-dpdf', '-painters');
    Handle_Figure1.Color = [1, 1, 1];
    
    Handle_Figure2.Color = 'none';
    set(Handle_Figure2, 'InvertHardCopy', 'off');
    print(Handle_Figure2, 'RealizedAW2_Saturation_Rejection_Disturbance', '-dpdf', '-painters');
    Handle_Figure2.Color = [1, 1, 1];
    
    % Handle_Figure3.Color = 'none';
    % set(Handle_Figure3, 'InvertHardCopy', 'off');
    % print(Handle_Figure3, 'RealizedAW2_Saturation_Rejection_AWVariables', '-dpdf', '-painters');
    % Handle_Figure3.Color = [1, 1, 1];
end
