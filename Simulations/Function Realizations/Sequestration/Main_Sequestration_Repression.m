%% Clear Workspace
close all
clear
clc

%% System 
% Stoichiometry Matrix and Propensity Function
StoichiometryMatrix = StoichiometryMatrix_Sequestration();
PropensityFunction = @PropensityFunction_Sequestration;
% System Parameters
Parameters.u_0 = 2;
Parameters.eta = 100;
Parameters.delta_1 = 1;
Parameters.delta_2 = 1;
Parameters.k = 1;
Parameters.kappa = 1;
% Input Parameters u(t) = u_max * (1 - exp(-a*t))
u_max_vector = [1, 3, 5, 7, 9];
a = 0.5;
% Sequestration
eta_vector = [1, 5, 20, 70];

%% Simulation Settings
tf = 15;               % Final Time
N_t = 1000;             % Time Samples
Solver = 'ODE23s';      % ODE Solver
IC = zeros(2,1);        % Initial Condition
                
%% Simulating the Full Model
t = cell(length(u_max_vector), 1);
V = cell(length(u_max_vector), 1);
for i = 1 : length(u_max_vector)
    Parameters.u = @(t) u_max_vector(i) * (1 - exp(-a*t));
    [t{i}, V{i}] = DSA_TimeVarying(StoichiometryMatrix, PropensityFunction, Parameters, IC, tf, N_t, Solver);
end

%% Static Input/Output Maps
f_eta = @(u, eta)( Parameters.k ./ (1 + (1/2/Parameters.kappa) * ((u - Parameters.u_0)/Parameters.delta_2 - Parameters.delta_1/eta + ...
                                      sqrt( ((u - Parameters.u_0)/Parameters.delta_2 - Parameters.delta_1/eta).^2 + ...
                                      4*Parameters.delta_1*u / Parameters.delta_2 / eta))) );
f_infinity = @(u)( Parameters.k ./ (1 + (1/Parameters.delta_2/Parameters.kappa) * max( (u - Parameters.u_0), 0 )) ); 

%% Generate Static Input/Output Maps
u_Grid = linspace(0, 10, 100);
y_eta = f_eta(u_Grid, Parameters.eta);
y_infinity = f_infinity(u_Grid);
y_bar_eta = f_eta(u_max_vector, Parameters.eta);
y_bar_infinity = f_infinity(u_max_vector);
y_eta_vector = zeros(length(eta_vector), length(u_Grid));
for i = 1 : length(eta_vector)
    y_eta_vector(i,:) = f_eta(u_Grid, eta_vector(i));
end

%% Figure Settings
ScalingFactor = 1;
SS = 4;
Figure_Width = 3 * SS;
Figure_Height = 3 * SS;
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

%% Set Figure 3
Handle_Figure3 = figure();
    Handle_Figure3.Color = [1 1 1];
    Handle_Figure3.PaperUnits = 'centimeters';
    Handle_Figure3.Units = 'centimeters';
    Handle_Figure3.Position = [0, 0, Figure_Width, Figure_Height];
    Handle_Figure3.PaperPositionMode = 'auto';
    Handle_Figure3.PaperSize = [Handle_Figure3.PaperPosition(3), Handle_Figure3.PaperPosition(4)];
    
%% Set Figure 4
Handle_Figure4 = figure();
    Handle_Figure4.Color = [1 1 1];
    Handle_Figure4.PaperUnits = 'centimeters';
    Handle_Figure4.Units = 'centimeters';
    Handle_Figure4.Position = [0, 0, Figure_Width, Figure_Height];
    Handle_Figure4.PaperPositionMode = 'auto';
    Handle_Figure4.PaperSize = [Handle_Figure4.PaperPosition(3), Handle_Figure4.PaperPosition(4)];

%% Axis for Output Response
Handle_Axis1 = axes(Handle_Figure1);
    Handle_Axis1.Box = 'on';
  	Handle_Axis1.BoxStyle = 'full';
 	Handle_Axis1.LineWidth = LineWidth_Thin;
  	Handle_Axis1.FontSize = FontSize;
 	hold(Handle_Axis1, 'on');
 	grid(Handle_Axis1, 'on');
   	Handle_Axis1.XMinorGrid = 'on';
  	Handle_Axis1.YMinorGrid = 'on';
   	Handle_Axis1.Title.String = 'Output Dynamic Response';
  	Handle_Axis1.XLabel.String = '$t$';
    Handle_Axis1.YLabel.String = '$y(t)$';
    Handle_Axis1.XLabel.Interpreter = 'latex';
    Handle_Axis1.YLabel.Interpreter = 'latex';
    Handle_Axis1.YLim = [0, y_eta(1)];
for i = 1 : length(u_max_vector)
    plot(Handle_Axis1, t{i}, Parameters.k ./ (1 + V{i}(2,:)/Parameters.kappa), 'Color', Colors(i,:), 'LineWidth', LineWidth_Thick);
end

%% Axis for Applied Input
Handle_Axis2 = axes(Handle_Figure2);
    Handle_Axis2.Box = 'on';
  	Handle_Axis2.BoxStyle = 'full';
 	Handle_Axis2.LineWidth = LineWidth_Thin;
  	Handle_Axis2.FontSize = FontSize;
 	hold(Handle_Axis2, 'on');
 	grid(Handle_Axis2, 'on');
   	Handle_Axis2.XMinorGrid = 'on';
  	Handle_Axis2.YMinorGrid = 'on';
   	Handle_Axis2.Title.String = 'Applied Input';
  	Handle_Axis2.XLabel.String = '$t$';
    Handle_Axis2.YLabel.String = '$u(t)$';
    Handle_Axis2.XLabel.Interpreter = 'latex';
    Handle_Axis2.YLabel.Interpreter = 'latex';
for i = 1 : length(u_max_vector)
    Parameters.u = @(t) u_max_vector(i) * (1 - exp(-a*t));
    plot(Handle_Axis2, t{i}, Parameters.u(t{i}), 'Color', Colors(i,:), 'LineWidth', LineWidth_Thick);
end
plot(Handle_Axis2, [0, tf], Parameters.u_0 * [1, 1], 'Color', 'k', 'LineWidth', LineWidth_Thin, 'LineStyle', '--');

%% Axis for Static Map
Handle_Axis3 = axes(Handle_Figure3);
    Handle_Axis3.Box = 'on';
  	Handle_Axis3.BoxStyle = 'full';
 	Handle_Axis3.LineWidth = LineWidth_Thin;
  	Handle_Axis3.FontSize = FontSize;
 	hold(Handle_Axis3, 'on');
 	grid(Handle_Axis3, 'on');
   	Handle_Axis3.XMinorGrid = 'on';
  	Handle_Axis3.YMinorGrid = 'on';
   	Handle_Axis3.Title.String = 'Input/Output Static Map';
    Handle_Axis3.XLabel.Interpreter = 'latex';
    Handle_Axis3.YLabel.Interpreter = 'latex';
  	Handle_Axis3.XLabel.String = '$\bar u$';
    Handle_Axis3.YLabel.String = '$\bar y$';
    Handle_Axis3.YLim = [0, y_eta(1)];
plot(Handle_Axis3, u_Grid, y_eta, 'Color', 'k', 'LineWidth', LineWidth_Thick);
for i = 1 : length(u_max_vector)
    plot(Handle_Axis3, u_max_vector(i), y_bar_infinity(i), 'MarkerFaceColor', Colors(i,:), 'MarkerSize', MarkerSize, 'Marker', 'o');
end

%% Axis for Convergence
Handle_Axis4 = axes(Handle_Figure4);
    Handle_Axis4.Position = [0.18, 0.16, 0.78, 0.75];
    Handle_Axis4.Box = 'on';
  	Handle_Axis4.BoxStyle = 'full';
 	Handle_Axis4.LineWidth = LineWidth_Thin;
  	Handle_Axis4.FontSize = FontSize;
 	hold(Handle_Axis4, 'on');
 	grid(Handle_Axis4, 'on');
   	Handle_Axis4.XMinorGrid = 'on';
  	Handle_Axis4.YMinorGrid = 'on';
   	Handle_Axis4.Title.String = 'Input/Output Static Map';
    Handle_Axis4.XLabel.Interpreter = 'latex';
    Handle_Axis4.YLabel.Interpreter = 'latex';
  	Handle_Axis4.XLabel.String = '$\bar u$';
    Handle_Axis4.YLabel.String = '$\bar y$';
    Handle_Axis4.YLim = [0, 1.1];
plot(Handle_Axis4, u_Grid, y_infinity, 'Color', 'k', 'LineWidth', LineWidth_Thick, 'LineStyle', '--');
for i = 1 : length(eta_vector)
    plot(Handle_Axis4, u_Grid, y_eta_vector(i,:), 'Color', 0.2 + 0.7 * [1, 1, 1] * (length(eta_vector) - i) / (length(eta_vector) - 1), 'LineWidth', LineWidth);
end

%% Save Figures
Save_Flag = 1;
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
    print(Handle_Figure1, 'Output_Seq_Rep', '-dpdf', '-painters');
    Handle_Figure1.Color = [1, 1, 1];
    
    Handle_Figure2.Color = 'none';
    set(Handle_Figure2, 'InvertHardCopy', 'off');
    print(Handle_Figure2, 'Input_Seq_Rep', '-dpdf', '-painters');
    Handle_Figure2.Color = [1, 1, 1];
    
    Handle_Figure3.Color = 'none';
    set(Handle_Figure3, 'InvertHardCopy', 'off');
    print(Handle_Figure3, 'Map_Seq_Rep', '-dpdf', '-painters');
    Handle_Figure3.Color = [1, 1, 1];
    
    Handle_Figure4.Color = 'none';
    set(Handle_Figure4, 'InvertHardCopy', 'off');
    print(Handle_Figure4, 'Convergence_Seq_Rep', '-dpdf', '-painters');
    Handle_Figure4.Color = [1, 1, 1];
end

