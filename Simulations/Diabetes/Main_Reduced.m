%% Clear Workspace
% close all
clear
clc
Save_Flag = 0;

%% Load Model
sbioloadproject('insulindemo', 'm1');
warnSettings = warning('off', 'SimBiology:DimAnalysisNotDone_MatlabFcn_Dimensionless');             % Suppress Warning
% addevent(m1, '[Plasma Glu] < 0*[Plasma Glu]', '[Plasma Glu] = 0 * [Plasma Glu]');

%% Simulation & Presentation Parameters
UB = 130;                               % [milligrams/decileter]    (Upper Bound)
LB = 80;                                % [milligrams/decileter]    (Lower Bound)
Setpoint = 100;                         % [milligrams/decileter]    (Lower Bound)
tf = 5036 + 48;                         % [hour]   
t0 = 5000;                              % [hour]
N_Cells = 1e9;                          % [dimensionless]

%% Disturbance
DisturbanceName = 'kp1';            
DisturbanceTime1 = 5024;                	% [hour]
DisturbanceFactor1 = 1.35/2.7;
DisturbanceTime2 = 5024 + 6;                	% [hour]
DisturbanceFactor2 = 1/DisturbanceFactor1;
addparameter(m1, 'DisturbanceTime1', DisturbanceTime1, 'Units', 'hour');
addparameter(m1, 'DisturbanceFactor1', DisturbanceFactor1, 'Units', 'dimensionless');
addparameter(m1, 'DisturbanceTime2', DisturbanceTime2, 'Units', 'hour');
addparameter(m1, 'DisturbanceFactor2', DisturbanceFactor2, 'Units', 'dimensionless');
Parameter = sbioselect(m1, 'Name', DisturbanceName);
    Parameter.Constant = false;
    Parameter.ConstantValue = false;
addevent(m1, 'time >= DisturbanceTime1', [DisturbanceName, ' = ', DisturbanceName, '*', num2str(DisturbanceFactor1)]);
addevent(m1, 'time >= DisturbanceTime2', [DisturbanceName, ' = ', DisturbanceName, '*', num2str(DisturbanceFactor2)]);

%% Saturation Parameters
kappa_a = 5e-9;  
kappa_s = 100;
kappa_P = 5e-9;

%% Integral Controller Parameters
mu = (100 / (6e23)) * 1e12 * 100;                	% [picomole/hour]
theta = (mu/Setpoint) * (1 + Setpoint/kappa_s);     % [(picomole/(milligram/deciliter))/hour]
eta = 1e10;                                         % [1/((milligram) * hour)]
k = 10 * N_Cells;                                   % [1/hour]
diffusion = (1/0.5);                                % [1/hour]
delta_c = 0;                                        % [1/hour]

addparameter(m1, 'mu', mu, 'Units', 'picomole/hour');
addparameter(m1, 'theta', theta, 'Units', '(picomole/(milligram/deciliter))/hour');
addparameter(m1, 'eta', eta, 'Units', '1/((picomole) * hour)');
addparameter(m1, 'k', k, 'Units', '1/hour');
addparameter(m1, 'diffusion', diffusion, 'Units', '1/hour');
addparameter(m1, 'delta_c', delta_c, 'Units', '1/hour');
addparameter(m1, 'Setpoint', Setpoint, 'Units', 'milligram/deciliter');
addparameter(m1, 'kappa_a', kappa_a, 'Units', 'picomole');
addparameter(m1, 'kappa_s', kappa_s, 'Units', 'milligram/deciliter');

%% Proportional Controller Parameters
alpha_Hill = mu/100;                                % [picomole/hour]
kappa_Hill = 200;                                   % [milligram/deciliter]
n_Hill = 4;                                         % [dimensionless]
delta = log(2)/0.5;                                 % [1/hour]
alpha_P = 0 * 700 * N_Cells;                        % [1/hour] 3*delta*50/mu 

addparameter(m1, 'kappa_P', kappa_P, 'Units', 'picomole');
addparameter(m1, 'alpha_Hill', alpha_Hill, 'Units', 'picomole/hour');
addparameter(m1, 'kappa_Hill', kappa_Hill, 'Units', 'milligram/deciliter');
addparameter(m1, 'n', n_Hill, 'Units', 'dimensionless');
addparameter(m1, 'delta', delta, 'Units', '1/hour');
addparameter(m1, 'alpha_P', alpha_P, 'Units', '1/hour');

%% Anti-Windup Parameters
alpha_1 = 1; beta_1 = 1e-8;
alpha_2 = 1; beta_2 = 0.1e-9;
addparameter(m1, 'alpha_1', alpha_1, 'Units', '1/hour');
addparameter(m1, 'alpha_2', alpha_2, 'Units', '1/hour');
addparameter(m1, 'beta_1', beta_1, 'Units', 'picomole');
addparameter(m1, 'beta_2', beta_2, 'Units', 'picomole');

%% Outputs to Plot
OutputNames = {'Z', 'Cell Ins'}; 

%% Glucose Dosage
mealDose = sbioselect(m1, 'Name', 'Single Meal');
    mealDose.Amount = 0;        % [gram]  
    mealDose.StartTime = t0;    % [hour]

%% Simulation Settings
configset = getconfigset(m1, 'active');
configset.StopTime = tf;
configset.SolverType = 'ode15s';
    
%% Normal Subject With No Control
% Simulation
NormalSim_NoControl = sbiosimulate(m1, configset, [], mealDose);                                   

%% Create Type-1 Diabetes Model
delete(m1.Reactions(10:11));

%% AIF Control Reactions
% Add Controller Compartment
Controller = addcompartment(m1, 'Controller', 'Units', 'liter');
% Controller Species
addspecies(Controller, 'Z', 0, 'Units', 'picomole');
addspecies(Controller, 'Z_3', 0, 'Units', 'picomole');
addspecies(Controller, 'Cell Ins', 0, 'Units', 'picomole');
% Reference Reaction
Reaction = addreaction(m1, 'null -> Controller.Z');
Reaction.ReactionRate = 'mu';
% Sensing Reaction
Reaction = addreaction(m1, 'Controller.Z -> null');
Reaction.ReactionRate = 'theta * [Glucose appearance].[Plasma Glu Conc] / (1 + [Glucose appearance].[Plasma Glu Conc]/kappa_s)';
% Transcription Reaction for P-Control
Reaction = addreaction(m1, 'null -> Controller.Z_3');
Reaction.ReactionRate = 'alpha_Hill * ([Glucose appearance].[Plasma Glu Conc] / kappa_Hill)^n / (([Glucose appearance].[Plasma Glu Conc] / kappa_Hill)^n + 1)';
% Degradation of P-Control Species
Reaction = addreaction(m1, 'Controller.Z_3 -> null');
Reaction.ReactionRate = 'delta * Controller.Z_3';
% Insulin Production Reaction in the Cell
Reaction = addreaction(m1, 'null -> Controller.[Cell Ins]');
Reaction.ReactionRate = 'k * max(-Controller.Z,0) / (1 + max(-Controller.Z, 0)/kappa_a) + alpha_P * Controller.Z_3 / (1 + Controller.Z_3/kappa_P)';
% Actuation Reaction
Reaction = addreaction(m1, 'Controller.[Cell Ins] -> [Insulin secretion].[Plasma Ins]');
Reaction.ReactionRate = 'diffusion * Controller.[Cell Ins]';
% Anti-Windup Reactions
Reaction = addreaction(m1, 'null -> Controller.Z');
Reaction.ReactionRate = 'alpha_1 * max(max(-Controller.Z, 0) - beta_1,0)';
Reaction = addreaction(m1, 'Controller.Z -> null');
Reaction.ReactionRate = 'alpha_2 * max(max(Controller.Z,0) - beta_2,0)';

%% AIF Control
% Simulation
DiabeticType1Sim_AW = sbiosimulate(m1, configset, [], mealDose);
param1 = sbioselect(m1, 'Type', 'parameter', 'Name', 'alpha_1');
param2 = sbioselect(m1, 'Type', 'parameter', 'Name', 'alpha_2');
param1.Value = 0;
param2.Value = 0;
DiabeticType1Sim = sbiosimulate(m1, configset, [], mealDose);

%% General Figure Settings
SS = 4; % Screen Scale
Figure_Width = 5 * SS;
Figure_Height = 5 * SS;
FontSize = 4 * SS;
FontSize_Small = 4 * SS;
LineWidth_Thin = 0.25 * SS;
LineWidth = 0.75 * SS;
MarkerSize = 3 * SS;

%% Color Choices
NormalColor = [0, 0, 0] / 255;              % Black
NCColor = [228,26,28] / 255;            	% Red
OLColor = [152,78,163] / 255;               % Purple
PIColor = [77,175,74] / 255;                % Green
PColor = [200, 200, 200] / 255;             % Gray
TubeColor = [77,175,74] / 255;          	% Green
IColor_AW = [55,126,184] / 255;             % Blue
IColor = [228,26,28] / 255;                 % Red

%% Figure 1 Settings
Handle_Figure1 = figure();
    Handle_Figure1.Color = [1 1 1];
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [-40, 40, 2.5*Figure_Width, Figure_Height];
    Handle_Figure1.PaperPositionMode = 'auto';
    Handle_Figure1.PaperSize = [Handle_Figure1.PaperPosition(3), Handle_Figure1.PaperPosition(4)];
figure(Handle_Figure1);

%% Figure 2 Settings
Handle_Figure2 = figure();
    Handle_Figure2.Color = [1 1 1];
    Handle_Figure2.PaperUnits = 'centimeters';
    Handle_Figure2.Units = 'centimeters';
    Handle_Figure2.Position = [-40, 40, Figure_Width, Figure_Height];
    Handle_Figure2.PaperPositionMode = 'auto';
    Handle_Figure2.PaperSize = [Handle_Figure2.PaperPosition(3), Handle_Figure2.PaperPosition(4)];
figure(Handle_Figure2);

%% Figure 1 Plots
set(0,'CurrentFigure', Handle_Figure1);
N_Rows = 3;
N_Columns = ceil(length(OutputNames) / N_Rows);
for i = 1:length(OutputNames)
    % Create Subplot
    Handle_Axis = subplot(N_Rows, N_Columns, i);
    hold(Handle_Axis, 'on');
    % Extract Simulations
    [tDiabetic_AW, yDiabetic_AW]  = DiabeticType1Sim_AW.selectbyname(OutputNames{i});
    [tDiabetic, yDiabetic]  = DiabeticType1Sim.selectbyname(OutputNames{i});
    plot(Handle_Axis, tDiabetic_AW, yDiabetic_AW, 'LineStyle', '-', 'LineWidth', LineWidth, 'Color', IColor_AW);
    plot(Handle_Axis, tDiabetic, yDiabetic, 'LineStyle', '-', 'LineWidth', LineWidth, 'Color', IColor);
    % Annotate Axis
    OutputParam = sbioselect(m1, 'Name', OutputNames{i});  
    Handle_Axis.Title.String = OutputNames{i};
    Handle_Axis.XLabel.String = 'Time (hr)';
    if strcmp(OutputParam.Type, 'parameter')
        Handle_Axis.YLabel.String = OutputParam.ValueUnits;
    else
        Handle_Axis.YLabel.String = OutputParam.InitialAmountUnits;
    end
    grid(Handle_Axis, 'on');
    Handle_Axis.Box = 'on';
	Handle_Axis.XLim = [t0-3, tf];
    Handle_Axis.XTick = linspace(t0, tf, 5);
    Handle_Axis.XTickLabel = linspace(0, tf-t0, 5);
end

%% Figure 2 Plots
t_Start = t0-3;
% Select Glucose
Output = 'Plasma Glu Conc';
% Create Axis
Handle_Axis2 = axes(Handle_Figure2);
	Handle_Axis2.Position = [0.1, 0.58, 0.86, 0.38];
	hold(Handle_Axis2, 'on');
	Handle_Axis2.Box = 'on';
 	Handle_Axis2.FontSize = FontSize;
  	grid(Handle_Axis2, 'on');
 	Handle_Axis2.XMinorGrid = 'off';
	Handle_Axis2.YMinorGrid = 'off';
 	Handle_Axis2.XLim = [t_Start, tf];
    Handle_Axis2.XTick = linspace(t0, tf, 7);
    Handle_Axis2.XTickLabel = linspace(0, tf-t0, 7);
    Handle_Axis2.YTick = [0, LB, UB, 200, 300];
    Handle_Axis2.Title.String = 'Plasma Glucose';
    Handle_Axis2.YLabel.String = 'mg/dl';
    Handle_Axis2.XLabel.String = 'Time (hr)';
%     Handle_Axis2.YLim = [0, 350];
% Extract Response of Glucose
[tNormal_NoControl, yNormal_NoControl]  = NormalSim_NoControl.selectbyname(Output);
[tDiabetic_AW, yDiabetic_AW]  = DiabeticType1Sim_AW.selectbyname(Output);
[tDiabetic, yDiabetic]  = DiabeticType1Sim.selectbyname(Output);
% Plot Responses
% plot(Handle_Axis2, [t_Start; tNormal_NoControl], [yNormal_NoControl(1); yNormal_NoControl], 'LineStyle', '-', 'LineWidth', LineWidth, 'Color', NormalColor);
plot(Handle_Axis2, [t_Start; tDiabetic_AW], [yDiabetic_AW(1); yDiabetic_AW], 'LineStyle', '-', 'LineWidth', LineWidth*0.8, 'Color', IColor_AW);
plot(Handle_Axis2, [t_Start; tDiabetic], [yDiabetic(1); yDiabetic], 'LineStyle', '-', 'LineWidth', LineWidth*0.8, 'Color', IColor);
plot(Handle_Axis2, [t_Start, tf], ((mu/theta) / (1 - mu/theta/kappa_s)) *[1, 1], 'LineStyle', '--', 'LineWidth', LineWidth, 'Color', 'k');
plot(Handle_Axis2, [t0, t0], Handle_Axis2.YLim, 'LineStyle', '--', 'LineWidth', LineWidth_Thin*2, 'Color', 0.5* [1, 1, 1]);
plot(Handle_Axis2, DisturbanceTime1*[1, 1], Handle_Axis2.YLim, 'LineStyle', '--', 'LineWidth', LineWidth_Thin, 'Color', [1, 0, 0]);
% Tube for Tolerance
plot(Handle_Axis2, [t_Start, tf], UB*[1, 1], 'LineStyle', '--', 'LineWidth', LineWidth_Thin, 'Color', 'k');
plot(Handle_Axis2, [t_Start, tf], LB*[1, 1], 'LineStyle', '--', 'LineWidth', LineWidth_Thin, 'Color', 'k');
Handle_Patch = patch(Handle_Axis2, [t_Start, tf, tf, t_Start], [LB, LB, UB, UB], 'green');
    Handle_Patch.FaceColor = TubeColor;
    Handle_Patch.EdgeColor = 'none';
    Handle_Patch.FaceAlpha = 0.2;
% Annotations
annotation(Handle_Figure2,'arrow',[0.166830026455025 0.166830026455025],...
    [0.971125220458554 0.961125220458554],...
    'Color',[0.650980392156863 0.650980392156863 0.650980392156863],...
    'LineWidth',10,...
    'HeadWidth',20,...
    'HeadStyle','plain',...
    'HeadLength',12);
annotation(Handle_Figure2,'textbox',...
    [0.162020061728395 0.916471119929453 0.0758472222222221 0.0493888888888886],...
    'String',{'Meal'},...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(Handle_Figure2,'arrow',[0.695930555555555 0.695930555555555],...
    [0.972888888888889 0.962888888888889],'Color',[1 0 0],'LineWidth',10,...
    'HeadWidth',20,...
    'HeadStyle','plain',...
    'HeadLength',12);
annotation(Handle_Figure2,'textbox',...
    [0.694669532627865 0.917110449735451 0.218910714285715 0.0487524250440903],...
    'String',{'Increase in EGP'},...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(Handle_Figure2,'textbox',...
    [0.754634259259258 0.89241909171076 0.0831082451499127 0.0487524250440903],...
    'String','Rate',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Select Insulin
Output = 'Plasma Ins Conc';
% Create Axis
Handle_Axis3 = axes(Handle_Figure2);
	Handle_Axis3.Position = [0.1, 0.08, 0.86, 0.38];
	hold(Handle_Axis3, 'on');
	Handle_Axis3.Box = 'on';
 	Handle_Axis3.FontSize = FontSize;
  	grid(Handle_Axis3, 'on');
 	Handle_Axis3.XMinorGrid = 'off';
	Handle_Axis3.YMinorGrid = 'off';
 	Handle_Axis3.XLim = [t_Start, tf];
    Handle_Axis3.XTick = linspace(t0, tf, 5);
    Handle_Axis3.XTickLabel = linspace(0, tf-t0, 5);
    Handle_Axis3.Title.String = 'Plasma Insulin';
    Handle_Axis3.XLabel.String = 'Time (hr)';
    Handle_Axis3.YLabel.String = 'pmol/l';
%     Handle_Axis3.YLim = [0, 125];
% Extract Response of Glucose
[tNormal_NoControl, yNormal_NoControl]  = NormalSim_NoControl.selectbyname(Output);
[tDiabetic_AW, yDiabetic_AW]  = DiabeticType1Sim_AW.selectbyname(Output);
[tDiabetic, yDiabetic]  = DiabeticType1Sim.selectbyname(Output);
% Plot Responses
% plot(Handle_Axis3, [t_Start; tNormal_NoControl], [yNormal_NoControl(1); yNormal_NoControl], 'LineStyle', '-', 'LineWidth', LineWidth, 'Color', NormalColor);
plot(Handle_Axis3, [t_Start; tDiabetic_AW], [yDiabetic_AW(1); yDiabetic_AW], 'LineStyle', '-', 'LineWidth', LineWidth*0.8, 'Color', IColor_AW);
plot(Handle_Axis3, [t_Start; tDiabetic], [yDiabetic(1); yDiabetic], 'LineStyle', '-', 'LineWidth', LineWidth*0.8, 'Color', IColor);
plot(Handle_Axis3, [t0, t0], Handle_Axis3.YLim, 'LineStyle', '--', 'LineWidth', LineWidth_Thin*2, 'Color', 0.5* [1, 1, 1]);
plot(Handle_Axis3, DisturbanceTime1*[1, 1], Handle_Axis3.YLim, 'LineStyle', '--', 'LineWidth', LineWidth_Thin, 'Color', [1, 0, 0]);

%% Save Figure
if Save_Flag == 1
    print(Handle_Figure2, 'Type1', '-dpdf');
end
