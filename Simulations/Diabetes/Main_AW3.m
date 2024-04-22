%% Clear Workspace
close all
clear
clc
Save_Flag = 1;

%% Load Model
sbioloadproject('insulindemo', 'm1');
warnSettings = warning('off', 'SimBiology:DimAnalysisNotDone_MatlabFcn_Dimensionless');             % Suppress Warning
% addevent(m1, '[Plasma Glu] < 0*[Plasma Glu]', '[Plasma Glu] = 0 * [Plasma Glu]');

%% Simulation & Presentation Parameters
UB = 126;                               % [milligrams/decileter]    (Upper Bound)
LB = 70;                                % [milligrams/decileter]    (Lower Bound)
Setpoint = 100;                         % [milligrams/decileter]    (Lower Bound)
tf = 5036 + 48;                         % [hour]   
t0 = 5000;                              % [hour]
N_Cells = 1e9;                          % [dimensionless]

%% Disturbance
DisturbanceName = 'kp1';            
DisturbanceTime1 = 5012;                	% [hour]
DisturbanceFactor1 = 1.35/2.7;
DisturbanceTime2 = 5012 + 6;                	% [hour]
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
alpha_P = 400 * N_Cells;                            % [1/hour] 3*delta*50/mu 

addparameter(m1, 'kappa_P', kappa_P, 'Units', 'picomole');
addparameter(m1, 'alpha_Hill', alpha_Hill, 'Units', 'picomole/hour');
addparameter(m1, 'kappa_Hill', kappa_Hill, 'Units', 'milligram/deciliter');
addparameter(m1, 'n', n_Hill, 'Units', 'dimensionless');
addparameter(m1, 'delta', delta, 'Units', '1/hour');
addparameter(m1, 'alpha_P', alpha_P, 'Units', '1/hour');

%% Anti-Windup Parameters
a_1 = 1; a_2 = 1;
delta_v1 = 1; delta_v2 = 1;
delta_w1 = 1; delta_w2 = 1;
kappa_v = 1e-8; kappa_w = 1e-8;
eta_v = 1e10;
eta_w = 1e10;
w_0 = 1e-8;
v_0 = 0.1e-9;

addparameter(m1, 'v_0', v_0, 'Units', 'picomole/hour');
addparameter(m1, 'w_0', w_0, 'Units', 'picomole/hour');
addparameter(m1, 'eta_v', eta_v, 'Units', '1/((picomole) * hour)');
addparameter(m1, 'eta_w', eta_w, 'Units', '1/((picomole) * hour)');
addparameter(m1, 'delta_v1', delta_v1, 'Units', '1/hour');
addparameter(m1, 'delta_v2', delta_v2, 'Units', '1/hour');
addparameter(m1, 'delta_w1', delta_w1, 'Units', '1/hour');
addparameter(m1, 'delta_w2', delta_w2, 'Units', '1/hour');
addparameter(m1, 'a_1', a_1, 'Units', '1/hour');
addparameter(m1, 'a_2', a_2, 'Units', '1/hour');
addparameter(m1, 'kappa_v', kappa_v, 'Units', 'picomole');
addparameter(m1, 'kappa_w', kappa_w, 'Units', 'picomole');

%% Outputs to Plot
OutputNames = {'Z_1', 'Z_2', 'Cell Ins'}; 

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
addspecies(Controller, 'Z_1', 0, 'Units', 'picomole');
addspecies(Controller, 'Z_2', 0, 'Units', 'picomole');
addspecies(Controller, 'Z_3', 0, 'Units', 'picomole');
addspecies(Controller, 'Cell Ins', 0, 'Units', 'picomole');
addspecies(Controller, 'V_1', 0, 'Units', 'picomole');
addspecies(Controller, 'V_2', 0, 'Units', 'picomole');
addspecies(Controller, 'W_1', 0, 'Units', 'picomole');
addspecies(Controller, 'W_2', 0, 'Units', 'picomole');
% Reference Reaction
Reaction = addreaction(m1, 'null -> Controller.Z_1');
Reaction.ReactionRate = 'mu * 1 / (1 + Controller.V_2/kappa_v)';
% Sensing Reaction
Reaction = addreaction(m1, 'null -> Controller.Z_2');
Reaction.ReactionRate = '(theta * [Glucose appearance].[Plasma Glu Conc] / (1 + [Glucose appearance].[Plasma Glu Conc]/kappa_s)) * (1 / (1 + Controller.W_2/kappa_w))';
% Sequestration Reaction
Reaction = addreaction(m1, 'Controller.Z_1 + Controller.Z_2 -> null');
Reaction.ReactionRate = 'eta * Controller.Z_1 * Controller.Z_2';
% Transcription Reaction for P-Control
Reaction = addreaction(m1, 'null -> Controller.Z_3');
Reaction.ReactionRate = 'alpha_Hill * ([Glucose appearance].[Plasma Glu Conc] / kappa_Hill)^n / (([Glucose appearance].[Plasma Glu Conc] / kappa_Hill)^n + 1)';
% Degradation of P-Control Species
Reaction = addreaction(m1, 'Controller.Z_3 -> null');
Reaction.ReactionRate = 'delta * Controller.Z_3';
% Insulin Production Reaction in the Cell
Reaction = addreaction(m1, 'null -> Controller.[Cell Ins]');
Reaction.ReactionRate = 'k * Controller.Z_2 / (1 + Controller.Z_2/kappa_a) + alpha_P * Controller.Z_3 / (1 + Controller.Z_3/kappa_P)';
% Actuation Reaction
Reaction = addreaction(m1, 'Controller.[Cell Ins] -> [Insulin secretion].[Plasma Ins]');
Reaction.ReactionRate = 'diffusion * Controller.[Cell Ins]';
% Dilution/Degradation Reactions
Reaction = addreaction(m1, 'Controller.Z_1 -> null');
Reaction.ReactionRate = 'delta_c * Controller.Z_1';
Reaction = addreaction(m1, 'Controller.Z_2 -> null');
Reaction.ReactionRate = 'delta_c * Controller.Z_2';
% Anti-Windup Reactions
Reaction = addreaction(m1, 'null -> Controller.V_1');
Reaction.ReactionRate = 'v_0';
Reaction = addreaction(m1, 'null -> Controller.V_2');
Reaction.ReactionRate = 'a_1 * Controller.Z_1';
Reaction = addreaction(m1, 'Controller.V_1 + Controller.V_2 -> null');
Reaction.ReactionRate = 'eta_v * Controller.V_1 * Controller.V_2';
Reaction = addreaction(m1, 'Controller.V_1 -> null');
Reaction.ReactionRate = 'delta_v1 * Controller.V_1';
Reaction = addreaction(m1, 'Controller.V_2 -> null');
Reaction.ReactionRate = 'delta_v2 * Controller.V_2';

Reaction = addreaction(m1, 'null -> Controller.W_1');
Reaction.ReactionRate = 'w_0';
Reaction = addreaction(m1, 'null -> Controller.W_2');
Reaction.ReactionRate = 'a_2 * Controller.Z_2';
Reaction = addreaction(m1, 'Controller.W_1 + Controller.W_2 -> null');
Reaction.ReactionRate = 'eta_w * Controller.W_1 * Controller.W_2';
Reaction = addreaction(m1, 'Controller.W_1 -> null');
Reaction.ReactionRate = 'delta_w1 * Controller.W_1';
Reaction = addreaction(m1, 'Controller.W_2 -> null');
Reaction.ReactionRate = 'delta_w2 * Controller.W_2';

%% AIF Control
% Simulation
DiabeticType1Sim_AW = sbiosimulate(m1, configset, [], mealDose);
param1 = sbioselect(m1, 'Type', 'parameter', 'Name', 'kappa_v');
param2 = sbioselect(m1, 'Type', 'parameter', 'Name', 'kappa_w');
param1.Value = inf;
param2.Value = inf;
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
NormalColor = [200, 200, 200] / 255;      	% Black
TubeColor = [77,175,74] / 255;          	% Green
IColor_AW = [55,126,184] / 255;             % Blue
IColor = [228,26,28] / 255;                 % Red
DisturbanceColor = [230, 0, 126] / 255;     % Pink
AW = [130, 55, 140] / 255;                  % Purple          

%% Figure 1 Settings
Handle_Figure1 = figure();
    Handle_Figure1.Color = [1 1 1];
    Handle_Figure1.PaperUnits = 'centimeters';
    Handle_Figure1.Units = 'centimeters';
    Handle_Figure1.Position = [-40, 40, Figure_Width, Figure_Height];
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
t_Start = t0-3;
% Select Z_1
Output = 'Z_1';
% Create Axis
Handle_Axis1 = axes(Handle_Figure1);
	Handle_Axis1.Position = [0.1, 0.58, 0.86, 0.38];
	hold(Handle_Axis1, 'on');
	Handle_Axis1.Box = 'on';
 	Handle_Axis1.FontSize = FontSize;
  	grid(Handle_Axis1, 'on');
 	Handle_Axis1.XMinorGrid = 'off';
	Handle_Axis1.YMinorGrid = 'off';
 	Handle_Axis1.XLim = [t_Start, tf];
    Handle_Axis1.XTick = linspace(t0, tf, 8);
    Handle_Axis1.XTickLabel = linspace(0, tf-t0, 8);
    Handle_Axis1.Title.String = 'Controller Species Z_1';
    Handle_Axis1.Title.Position(2) = 5.9e-8;
    Handle_Axis1.YLabel.String = 'pmol/cell';
    Handle_Axis1.XLabel.String = 'Time (hr)';
    Handle_Axis1.YLim = [0, 0.6e-7];
% Extract Response of Z_1
[tDiabetic_AW, yDiabetic_AW]  = DiabeticType1Sim_AW.selectbyname(Output);
[tDiabetic, yDiabetic]  = DiabeticType1Sim.selectbyname(Output);
% Plot Responses
plot(Handle_Axis1, [t_Start; tDiabetic_AW], [yDiabetic_AW(1); yDiabetic_AW], 'LineStyle', '-', 'LineWidth', LineWidth*0.8, 'Color', IColor_AW);
plot(Handle_Axis1, [t_Start; tDiabetic], [yDiabetic(1); yDiabetic], 'LineStyle', '-', 'LineWidth', LineWidth*0.8, 'Color', IColor);
plot(Handle_Axis1, DisturbanceTime1*[1, 1], Handle_Axis1.YLim, 'LineStyle', '-', 'LineWidth', LineWidth_Thin, 'Color', [1, 0, 0]);
plot(Handle_Axis1, DisturbanceTime2*[1, 1], Handle_Axis1.YLim, 'LineStyle', '-', 'LineWidth', LineWidth_Thin, 'Color', [0, 0, 1]);
% Tube for Disturbance
Handle_Patch = patch(Handle_Axis1, [DisturbanceTime1, DisturbanceTime2, DisturbanceTime2, DisturbanceTime1], [Handle_Axis1.YLim(1), Handle_Axis1.YLim(1), Handle_Axis1.YLim(2), Handle_Axis1.YLim(2)], DisturbanceColor);
    Handle_Patch.FaceColor = DisturbanceColor;
    Handle_Patch.EdgeColor = 'none';
    Handle_Patch.FaceAlpha = 0.2;
    
% Select V_2
Output = 'V_2';
% Create Axis
Handle_Axis11 = axes(Handle_Figure1);
	Handle_Axis11.Position = [0.5, 0.79, 0.435, 0.11];
	hold(Handle_Axis11, 'on');
	Handle_Axis11.Box = 'on';
 	Handle_Axis11.FontSize = FontSize/1.1;
  	grid(Handle_Axis11, 'on');
 	Handle_Axis11.XMinorGrid = 'off';
	Handle_Axis11.YMinorGrid = 'off';
 	Handle_Axis11.XLim = [t_Start, tf];
    Handle_Axis11.XTick = linspace(t0, tf, 8);
    Handle_Axis11.XTickLabel = linspace(0, tf-t0, 8);
    Handle_Axis11.Title.String = 'Anti-Windup Species V_2';
    Handle_Axis11.Title.Position(1:2) = [5.045e3, 1.15e-8];
%     Handle_Axis11.YLabel.String = 'pmol/cell';
    Handle_Axis11.XLabel.String = 'Time (hr)';
    Handle_Axis11.YLim = [0, 1.2e-8];
% Extract Response of W_2
[tDiabetic_AW, yDiabetic_AW]  = DiabeticType1Sim_AW.selectbyname(Output);
% Plot Responses
plot(Handle_Axis11, [t_Start; tDiabetic_AW], [yDiabetic_AW(1); yDiabetic_AW], 'LineStyle', '-', 'LineWidth', LineWidth*0.8, 'Color', AW);
plot(Handle_Axis11, DisturbanceTime1*[1, 1], Handle_Axis1.YLim, 'LineStyle', '-', 'LineWidth', LineWidth_Thin, 'Color', [1, 0, 0]);
plot(Handle_Axis11, DisturbanceTime2*[1, 1], Handle_Axis1.YLim, 'LineStyle', '-', 'LineWidth', LineWidth_Thin, 'Color', [0, 0, 1]);
% % Tube for Disturbance
Handle_Patch = patch(Handle_Axis11, [DisturbanceTime1, DisturbanceTime2, DisturbanceTime2, DisturbanceTime1], [Handle_Axis11.YLim(1), Handle_Axis11.YLim(1), Handle_Axis11.YLim(2), Handle_Axis11.YLim(2)], DisturbanceColor);
    Handle_Patch.FaceColor = DisturbanceColor;
    Handle_Patch.EdgeColor = 'none';
    Handle_Patch.FaceAlpha = 0.2;
    
% Select Z_2
Output = 'Z_2';
% Create Axis
Handle_Axis12 = axes(Handle_Figure1);
	Handle_Axis12.Position = [0.1, 0.08, 0.86, 0.38];
	hold(Handle_Axis12, 'on');
	Handle_Axis12.Box = 'on';
 	Handle_Axis12.FontSize = FontSize;
  	grid(Handle_Axis12, 'on');
 	Handle_Axis12.XMinorGrid = 'off';
	Handle_Axis12.YMinorGrid = 'off';
 	Handle_Axis12.XLim = [t_Start, tf];
    Handle_Axis12.XTick = linspace(t0, tf, 8);
    Handle_Axis12.XTickLabel = linspace(0, tf-t0, 8);
    Handle_Axis12.Title.String = 'Controller Species Z_2';
    Handle_Axis12.Title.Position(2) = 1.96e-8;
    Handle_Axis12.YLabel.String = 'pmol/cell';
    Handle_Axis12.XLabel.String = 'Time (hr)';
    Handle_Axis12.YLim = [0, 0.2e-7];
% Extract Response of Z_2
[tDiabetic_AW, yDiabetic_AW]  = DiabeticType1Sim_AW.selectbyname(Output);
[tDiabetic, yDiabetic]  = DiabeticType1Sim.selectbyname(Output);
% Plot Responses
plot(Handle_Axis12, [t_Start; tDiabetic_AW], [yDiabetic_AW(1); yDiabetic_AW], 'LineStyle', '-', 'LineWidth', LineWidth*0.8, 'Color', IColor_AW);
plot(Handle_Axis12, [t_Start; tDiabetic], [yDiabetic(1); yDiabetic], 'LineStyle', '-', 'LineWidth', LineWidth*0.8, 'Color', IColor);
plot(Handle_Axis12, DisturbanceTime1*[1, 1], Handle_Axis12.YLim, 'LineStyle', '-', 'LineWidth', LineWidth_Thin, 'Color', [1, 0, 0]);
plot(Handle_Axis12, DisturbanceTime2*[1, 1], Handle_Axis12.YLim, 'LineStyle', '-', 'LineWidth', LineWidth_Thin, 'Color', [0, 0, 1]);
% Tube for Disturbance
Handle_Patch = patch(Handle_Axis12, [DisturbanceTime1, DisturbanceTime2, DisturbanceTime2, DisturbanceTime1], [Handle_Axis12.YLim(1), Handle_Axis12.YLim(1), Handle_Axis12.YLim(2), Handle_Axis12.YLim(2)], DisturbanceColor);
    Handle_Patch.FaceColor = DisturbanceColor;
    Handle_Patch.EdgeColor = 'none';
    Handle_Patch.FaceAlpha = 0.2;

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
    Handle_Axis2.XTick = linspace(t0, tf, 8);
    Handle_Axis2.XTickLabel = linspace(0, tf-t0, 8);
    Handle_Axis2.YTick = [0, LB, UB, 200, 300];
    Handle_Axis2.Title.String = 'Plasma Glucose';
    Handle_Axis2.YLabel.String = 'mg/dl';
    Handle_Axis2.XLabel.String = 'Time (hr)';
    Handle_Axis2.YLim = [0, 200];
% Extract Response of Glucose
[tNormal_NoControl, yNormal_NoControl]  = NormalSim_NoControl.selectbyname(Output);
[tDiabetic_AW, yDiabetic_AW]  = DiabeticType1Sim_AW.selectbyname(Output);
[tDiabetic, yDiabetic]  = DiabeticType1Sim.selectbyname(Output);
% Plot Responses
Handle_Normal = plot(Handle_Axis2, [t_Start; tNormal_NoControl], [yNormal_NoControl(1); yNormal_NoControl], 'LineStyle', '-', 'LineWidth', LineWidth, 'Color', NormalColor);
Handle_AW = plot(Handle_Axis2, [t_Start; tDiabetic_AW], [yDiabetic_AW(1); yDiabetic_AW], 'LineStyle', '-', 'LineWidth', LineWidth*0.8, 'Color', IColor_AW);
Handle_NoAW = plot(Handle_Axis2, [t_Start; tDiabetic], [yDiabetic(1); yDiabetic], 'LineStyle', '-', 'LineWidth', LineWidth*0.8, 'Color', IColor);
Handle_Ref = plot(Handle_Axis2, [t_Start, tf], ((mu/theta) / (1 - mu/theta/kappa_s)) *[1, 1], 'LineStyle', '--', 'LineWidth', LineWidth, 'Color', 'k');
plot(Handle_Axis2, DisturbanceTime1*[1, 1], Handle_Axis2.YLim, 'LineStyle', '-', 'LineWidth', LineWidth_Thin, 'Color', [1, 0, 0]);
plot(Handle_Axis2, DisturbanceTime2*[1, 1], Handle_Axis2.YLim, 'LineStyle', '-', 'LineWidth', LineWidth_Thin, 'Color', [0, 0, 1]);
% Tube for Tolerance
plot(Handle_Axis2, [t_Start, tf], UB*[1, 1], 'LineStyle', '--', 'LineWidth', LineWidth_Thin, 'Color', 'k');
plot(Handle_Axis2, [t_Start, tf], LB*[1, 1], 'LineStyle', '--', 'LineWidth', LineWidth_Thin, 'Color', 'k');
Handle_Patch = patch(Handle_Axis2, [t_Start, tf, tf, t_Start], [LB, LB, UB, UB], 'green');
    Handle_Patch.FaceColor = TubeColor;
    Handle_Patch.EdgeColor = 'none';
    Handle_Patch.FaceAlpha = 0.2;
% Tube for Disturbance
Handle_Patch = patch(Handle_Axis2, [DisturbanceTime1, DisturbanceTime2, DisturbanceTime2, DisturbanceTime1], [Handle_Axis2.YLim(1), Handle_Axis2.YLim(1), Handle_Axis2.YLim(2), Handle_Axis2.YLim(2)], DisturbanceColor);
    Handle_Patch.FaceColor = DisturbanceColor;
    Handle_Patch.EdgeColor = 'none';
    Handle_Patch.FaceAlpha = 0.2;
% Legend
Handle_Legend = legend(Handle_Axis2, [Handle_NoAW, Handle_AW, Handle_Ref, Handle_Normal], ...
                                      {'Diabetic - Without Anti-Windup', 'Diabetic - With Anti-Windup', 'Setpoint', 'Normal'});
    Handle_Legend.Location = 'south east';
% Annotations
annotation(Handle_Figure2,'arrow',[0.247958774250439 0.247958774250439],...
    [0.971125220458554 0.961125220458554],'Color',[1 0 0],'LineWidth',10,...
    'HeadWidth',20,...
    'HeadStyle','plain',...
    'HeadLength',12);
annotation(Handle_Figure2,'textbox',...
    [0.0597272927689594 0.838875881834213 0.205026455026455 0.112874779541446],...
    'Color',[1 0 0],...
    'String',{'Extreme Drop','in EGP Rate','(50%)'},...
    'HorizontalAlignment','center',...
    'FontSize',15,...
    'EdgeColor','none');
annotation(Handle_Figure2,'arrow',[0.307923500881833 0.307923500881833],...
    [0.971125220458554 0.961125220458554],'Color',[0 0 1],'LineWidth',10,...
    'HeadWidth',20,...
    'HeadStyle','plain',...
    'HeadLength',12);
annotation(Handle_Figure2,'textbox',...
    [0.297844135802468 0.901237433862434 0.218910714285715 0.0487524250440904],...
    'Color',[0 0 1],...
    'String',{'EGP Rate Back','to Normal'},...
    'HorizontalAlignment','center',...
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
    Handle_Axis3.XTick = linspace(t0, tf, 8);
    Handle_Axis3.XTickLabel = linspace(0, tf-t0, 8);
    Handle_Axis3.Title.String = 'Plasma Insulin';
    Handle_Axis3.XLabel.String = 'Time (hr)';
    Handle_Axis3.YLabel.String = 'pmol/l';
    Handle_Axis3.YLim = [0, 50];
% Extract Response of Glucose
[tNormal_NoControl, yNormal_NoControl]  = NormalSim_NoControl.selectbyname(Output);
[tDiabetic_AW, yDiabetic_AW]  = DiabeticType1Sim_AW.selectbyname(Output);
[tDiabetic, yDiabetic]  = DiabeticType1Sim.selectbyname(Output);
% Plot Responses
plot(Handle_Axis3, [t_Start; tNormal_NoControl], [yNormal_NoControl(1); yNormal_NoControl], 'LineStyle', '-', 'LineWidth', LineWidth, 'Color', NormalColor);
plot(Handle_Axis3, [t_Start; tDiabetic_AW], [yDiabetic_AW(1); yDiabetic_AW], 'LineStyle', '-', 'LineWidth', LineWidth*0.8, 'Color', IColor_AW);
plot(Handle_Axis3, [t_Start; tDiabetic], [yDiabetic(1); yDiabetic], 'LineStyle', '-', 'LineWidth', LineWidth*0.8, 'Color', IColor);
plot(Handle_Axis3, DisturbanceTime1*[1, 1], Handle_Axis3.YLim, 'LineStyle', '-', 'LineWidth', LineWidth_Thin, 'Color', [1, 0, 0]);
plot(Handle_Axis3, DisturbanceTime2*[1, 1], Handle_Axis3.YLim, 'LineStyle', '-', 'LineWidth', LineWidth_Thin, 'Color', [0, 0, 1]);
Handle_Patch = patch(Handle_Axis3, [DisturbanceTime1, DisturbanceTime2, DisturbanceTime2, DisturbanceTime1], [Handle_Axis2.YLim(1), Handle_Axis2.YLim(1), Handle_Axis2.YLim(2), Handle_Axis2.YLim(2)], DisturbanceColor);
    Handle_Patch.FaceColor = DisturbanceColor;
    Handle_Patch.EdgeColor = 'none';
    Handle_Patch.FaceAlpha = 0.2;
    
%% Save Figure
if Save_Flag == 1
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
    print(Handle_Figure1, 'Diabetes_Control_AW3', '-dpdf', '-painters');
    Handle_Figure1.Color = [1, 1, 1];
    
    Handle_Figure2.Color = 'none';
    set(Handle_Figure2, 'InvertHardCopy', 'off');
    print(Handle_Figure2, 'Diabetes_Response_AW3', '-dpdf', '-painters');
    Handle_Figure2.Color = [1, 1, 1];
end
