%% Clear Workspace
close all
clear
clc

%% Load
load GeneExp_Deterministic_Sweep_100
clearvars -except Error gamma_1_vector gamma_2_vector eta_vector

%% Truncate due to Numerical Errors at Large eta
Error(:,:,end-1:end) = [];
eta_vector(end-1:end) = [];

%% Simulations for Nominal Values
tf = 300;
t_Disturbance = [100, 150];
DisturbanceFactor = [2, 1/3];
L = 2;
% Nominal Values 1
gamma_1_Nominal1 = 1;
gamma_2_Nominal1 = 1;
t1 = cell(1, length(eta_vector));
X1 = cell(1, length(eta_vector));
Error_Nominal1 = zeros(1, length(eta_vector));
for i = 1 : length(eta_vector)
    [t1{i}, X1{i}, X_Reduced_Mapped1, Setpoint1, Error_Nominal1(i)] = GeneExp_Deterministic(eta_vector(i), gamma_1_Nominal1, gamma_2_Nominal1, tf, t_Disturbance, DisturbanceFactor);
end

% Nominal Values 2
gamma_1_Nominal2 = 0.5;
gamma_2_Nominal2 = 0.5;
t2 = cell(1, length(eta_vector));
X2 = cell(1, length(eta_vector));
Error_Nominal2 = zeros(1, length(eta_vector));
for i = 1 : length(eta_vector)
    [t2{i}, X2{i}, X_Reduced_Mapped2, Setpoint2, Error_Nominal2(i)] = GeneExp_Deterministic(eta_vector(i), gamma_1_Nominal2, gamma_2_Nominal2, tf, t_Disturbance, DisturbanceFactor);
end

% Nominal Values 3
gamma_1_Nominal3 = 3;
gamma_2_Nominal3 = 3;
t3 = cell(1, length(eta_vector));
X3 = cell(1, length(eta_vector));
Error_Nominal3 = zeros(1, length(eta_vector));
for i = 1 : length(eta_vector)
    [t3{i}, X3{i}, X_Reduced_Mapped3, Setpoint3, Error_Nominal3(i)] = GeneExp_Deterministic(eta_vector(i), gamma_1_Nominal3, gamma_2_Nominal3, tf, t_Disturbance, DisturbanceFactor);
end

%% Reshaping
n1 = size(Error, 1);
n2 = size(Error, 2);
n3 = size(Error, 3);
Error = permute(Error, [1, 3, 2]);
Error_Summary = 0 * eta_vector;
for i = 1 : length(eta_vector)
    Error_Summary(i) = mean(nanmean(squeeze(Error(:,i,:)))); 
end

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
Opacity = 0.25;
Colors = [39, 52, 139; ... % Blue
          227, 6, 19; ... % Red
          0, 152, 58; ... % Green
          128, 0, 128; ... % Purple
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
      
%% 3D Instensity Plot   
set(0, 'CurrentFigure', Handle_Figure1);
[x, y, z] = meshgrid(eta_vector, gamma_1_vector, gamma_2_vector);
Handle_Axis1 = subplot(2,1,1); 
hold(Handle_Axis1, 'on');
Handle_Slice = slice(Handle_Axis1, x, y, z, log10(Error), eta_vector, [], []);
for k = 1 : length(eta_vector)
    Handle_Slice(k).FaceAlpha = 0.6;
end
shading interp; 
% Colormap
colormap(Handle_Axis1, 'jet');
Handle_Colorbar = colorbar;
Handle_Colorbar.Ticks = log10([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]); 
Handle_Colorbar.TickLabels = {'10^{-6}', '10^{-5}', '10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}', '1'}; 
% Axis Attributes
Handle_Axis1.XScale = 'log';
Handle_Axis1.XLabel.String = '$\eta$';
grid(Handle_Axis1, 'on');
Handle_Axis1.YLabel.String = '$\gamma_1$';
Handle_Axis1.ZLabel.String = '$\gamma_2$';
Handle_Axis1.XLabel.Interpreter = 'latex';
Handle_Axis1.YLabel.Interpreter = 'latex';
Handle_Axis1.ZLabel.Interpreter = 'latex';
Handle_Axis1.View = [6.8, 25.5];
Handle_Axis1.FontSize = FontSize;
Handle_Axis1.XLim = [eta_vector(1), eta_vector(end)];
Handle_Axis1.YLim = [gamma_1_vector(1), gamma_1_vector(end)];
Handle_Axis1.ZLim = [gamma_2_vector(1), gamma_2_vector(end)];
Handle_Axis1.Title.String = 'Relative Error';
Handle_Axis1.XTick = eta_vector;
Handle_Axis1.YTick = linspace(1, 5, 5);
Handle_Axis1.ZTick = linspace(1, 5, 5);
Handle_Axis1.XLabel.Position = [12.0000, 0.1937, 0.2500];
Handle_Axis1.YLabel.Position(3) = 0.45;
% Nominal Points
for k = 1 : length(eta_vector)
    scatter3(Handle_Axis1, eta_vector(k), gamma_1_Nominal1, gamma_2_Nominal1, 'filled', 'MarkerFaceAlpha', Opacity + (1 - Opacity) * (k-1)/(length(eta_vector) - 1), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', Colors(1,:), 'SizeData', 400);
    scatter3(Handle_Axis1, eta_vector(k), gamma_1_Nominal2, gamma_2_Nominal2, 'filled', 'MarkerFaceAlpha', Opacity + (1 - Opacity) * (k-1)/(length(eta_vector) - 1), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', Colors(3,:), 'SizeData', 400);
    scatter3(Handle_Axis1, eta_vector(k), gamma_1_Nominal3, gamma_2_Nominal3, 'filled', 'MarkerFaceAlpha', Opacity + (1 - Opacity) * (k-1)/(length(eta_vector) - 1), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', Colors(4,:), 'SizeData', 400);
end

%% Error Summary
Handle_Axis2 = subplot(2,1,2);
hold(Handle_Axis2, 'on');
plot(Handle_Axis2, eta_vector, Error_Summary, 'LineWidth', LineWidth, 'Color', 'k');
plot(Handle_Axis2, eta_vector, Error_Nominal1, 'LineWidth', LineWidth, 'Color', Colors(1,:));
plot(Handle_Axis2, eta_vector, Error_Nominal2, 'LineWidth', LineWidth, 'Color', Colors(3,:));
plot(Handle_Axis2, eta_vector, Error_Nominal3, 'LineWidth', LineWidth, 'Color', Colors(4,:));
Handle_Axis2.XScale = 'log';
Handle_Axis2.FontSize = FontSize;
Handle_Axis2.XLabel.String = '$\eta$';
Handle_Axis2.XLabel.Interpreter = 'latex';
Handle_Axis2.YLabel.String = 'Relative Error';
grid(Handle_Axis2, 'on');
Handle_Legend = legend(Handle_Axis2, {'Mean', ...
                                      ['$(\gamma_1, \gamma_2) = (', num2str(gamma_1_Nominal1), ',', num2str(gamma_2_Nominal1), ')$'], ...
                                      ['$(\gamma_1, \gamma_2) = (', num2str(gamma_1_Nominal2), ',', num2str(gamma_2_Nominal2), ')$'], ...
                                      ['$(\gamma_1, \gamma_2) = (', num2str(gamma_1_Nominal3), ',', num2str(gamma_2_Nominal3), ')$'], ...
                                      });
Handle_Legend.Interpreter = 'latex';
                               
%% Nominal Simulations 1
set(0, 'CurrentFigure', Handle_Figure2);
Plot_Names = {'x_1', 'x_2', 'z_1', 'z_2'};
rows = 2;
columns = 2;
Handle_Axis3 = gobjects(rows, columns);
for i = 1 : rows
    for j = 1 : columns
        Handle_Axis3(i,j) = subplot(rows, columns, (i - 1)*columns + j); 
        Handle_Axis3(i,j).Box = 'on';
        Handle_Axis3(i,j).BoxStyle = 'full';
        Handle_Axis3(i,j).LineWidth = LineWidth_Thin;
        Handle_Axis3(i,j).FontSize = FontSize;
        hold(Handle_Axis3(i,j), 'on');
        grid(Handle_Axis3(i,j), 'on');
        Handle_Axis3(i,j).XMinorGrid = 'on';
        Handle_Axis3(i,j).YMinorGrid = 'on';
        Handle_Axis3(i,j).XLabel.String = '$t$';
        Handle_Axis3(i,j).XLabel.Interpreter = 'latex';
        Handle_Axis3(i,j).YLabel.String = ['$', Plot_Names{(i - 1)*columns + j}, '(t)$'];
        Handle_Axis3(i,j).YLabel.Interpreter = 'latex';
        for k = 1 : length(eta_vector)
            plot(Handle_Axis3(i,j), t1{k}, X1{k}((i - 1)*columns + j,:), 'Color', [Colors(1,:), Opacity + (1 - Opacity) * (k-1)/(length(eta_vector) - 1)], 'LineWidth', LineWidth_Thick/2, 'LineStyle', '-');
        end
        plot(Handle_Axis3(i,j), t1{k}, X_Reduced_Mapped1((i - 1)*columns + j,:), 'Color', Colors(2,:), 'LineWidth', LineWidth_Thick/2, 'LineStyle', '--');
        if (i - 1)*columns + j == L
            Handle_Setpoint1 = plot(Handle_Axis3(i,j), [0, t_Disturbance(1), t_Disturbance(1), t_Disturbance(2), t_Disturbance(2), tf], ...
                                                       [Setpoint1(1), Setpoint1(1), Setpoint1(2), Setpoint1(2), Setpoint1(3), Setpoint1(3)], ...
                                                       'LineWidth', LineWidth, 'Color', 'k', 'LineStyle', '--');
            legend(Handle_Axis3(i,j), Handle_Setpoint1, 'Setpoint');
        end
    end
end
Handle_Legend1 = legend(Handle_Axis3(rows, columns), [repmat('$\eta = $', length(eta_vector), 1), num2str(eta_vector')]);
Handle_Legend1 = legend(Handle_Axis3(rows, columns), [Handle_Legend1.String, 'Reduced Model']);
Handle_Legend1.Interpreter = 'latex';
Handle_Legend1.Location = 'Best';

%% Nominal Simulations 2
set(0, 'CurrentFigure', Handle_Figure3);
Plot_Names = {'x_1', 'x_2', 'z_1', 'z_2'};
rows = 2;
columns = 2;
Handle_Axis4 = gobjects(rows, columns);
for i = 1 : rows
    for j = 1 : columns
        Handle_Axis4(i,j) = subplot(rows, columns, (i - 1)*columns + j); 
        Handle_Axis4(i,j).Box = 'on';
        Handle_Axis4(i,j).BoxStyle = 'full';
        Handle_Axis4(i,j).LineWidth = LineWidth_Thin;
        Handle_Axis4(i,j).FontSize = FontSize;
        hold(Handle_Axis4(i,j), 'on');
        grid(Handle_Axis4(i,j), 'on');
        Handle_Axis4(i,j).XMinorGrid = 'on';
        Handle_Axis4(i,j).YMinorGrid = 'on';
        Handle_Axis4(i,j).XLabel.String = '$t$';
        Handle_Axis4(i,j).XLabel.Interpreter = 'latex';
        Handle_Axis4(i,j).YLabel.String = ['$', Plot_Names{(i - 1)*columns + j}, '(t)$'];
        Handle_Axis4(i,j).YLabel.Interpreter = 'latex';
        for k = 1 : length(eta_vector)
            plot(Handle_Axis4(i,j), t2{k}, X2{k}((i - 1)*columns + j,:), 'Color', [Colors(3,:), Opacity + (1 - Opacity) * (k-1)/(length(eta_vector) - 1)], 'LineWidth', LineWidth_Thick/2, 'LineStyle', '-');
        end
        plot(Handle_Axis4(i,j), t2{k}, X_Reduced_Mapped2((i - 1)*columns + j,:), 'Color', Colors(2,:), 'LineWidth', LineWidth_Thick/2, 'LineStyle', '--');
        if (i - 1)*columns + j == L
            Handle_Setpoint2 = plot(Handle_Axis4(i,j), [0, t_Disturbance(1), t_Disturbance(1), t_Disturbance(2), t_Disturbance(2), tf], ...
                                                       [Setpoint1(1), Setpoint1(1), Setpoint1(2), Setpoint1(2), Setpoint1(3), Setpoint1(3)], ...
                                                       'LineWidth', LineWidth, 'Color', 'k', 'LineStyle', '--');
            legend(Handle_Axis4(i,j), Handle_Setpoint2, 'Setpoint');
        end
    end
end
Handle_Legend2 = legend(Handle_Axis4(rows, columns), [repmat('$\eta = $', length(eta_vector), 1), num2str(eta_vector')]);
Handle_Legend2 = legend(Handle_Axis4(rows, columns), [Handle_Legend2.String, 'Reduced Model']);
Handle_Legend2.Interpreter = 'latex';
Handle_Legend2.Location = 'Best';

%% Nominal Simulations 3
set(0, 'CurrentFigure', Handle_Figure4);
Plot_Names = {'x_1', 'x_2', 'z_1', 'z_2'};
rows = 2;
columns = 2;
Handle_Axis5 = gobjects(rows, columns);
for i = 1 : rows
    for j = 1 : columns
        Handle_Axis5(i,j) = subplot(rows, columns, (i - 1)*columns + j); 
        Handle_Axis5(i,j).Box = 'on';
        Handle_Axis5(i,j).BoxStyle = 'full';
        Handle_Axis5(i,j).LineWidth = LineWidth_Thin;
        Handle_Axis5(i,j).FontSize = FontSize;
        hold(Handle_Axis5(i,j), 'on');
        grid(Handle_Axis5(i,j), 'on');
        Handle_Axis5(i,j).XMinorGrid = 'on';
        Handle_Axis5(i,j).YMinorGrid = 'on';
        Handle_Axis5(i,j).XLabel.String = '$t$';
        Handle_Axis5(i,j).XLabel.Interpreter = 'latex';
        Handle_Axis5(i,j).YLabel.String = ['$', Plot_Names{(i - 1)*columns + j}, '(t)$'];
        Handle_Axis5(i,j).YLabel.Interpreter = 'latex';
        for k = 1 : length(eta_vector)
            plot(Handle_Axis5(i,j), t3{k}, X3{k}((i - 1)*columns + j,:), 'Color', [Colors(4,:), Opacity + (1 - Opacity) * (k-1)/(length(eta_vector) - 1)], 'LineWidth', LineWidth_Thick/2, 'LineStyle', '-');
        end
        plot(Handle_Axis5(i,j), t3{k}, X_Reduced_Mapped3((i - 1)*columns + j,:), 'Color', Colors(2,:), 'LineWidth', LineWidth_Thick/2, 'LineStyle', '--');
        if (i - 1)*columns + j == L
            Handle_Setpoint3 = plot(Handle_Axis5(i,j), [0, t_Disturbance(1), t_Disturbance(1), t_Disturbance(2), t_Disturbance(2), tf], ...
                                                       [Setpoint1(1), Setpoint1(1), Setpoint1(2), Setpoint1(2), Setpoint1(3), Setpoint1(3)], ...
                                                       'LineWidth', LineWidth, 'Color', 'k', 'LineStyle', '--');
            legend(Handle_Axis5(i,j), Handle_Setpoint3, 'Setpoint');
        end
    end
end
Handle_Legend3 = legend(Handle_Axis5(rows, columns), [repmat('$\eta = $', length(eta_vector), 1), num2str(eta_vector')]);
Handle_Legend3 = legend(Handle_Axis5(rows, columns), [Handle_Legend3.String, 'Reduced Model']);
Handle_Legend3.Interpreter = 'latex';
Handle_Legend3.Location = 'Best';

%% Positioning of the Axes
Handle_Axis1.Position = [0.05, 0.45, 0.85, 0.5];
Handle_Colorbar.Position(1) = 0.95;
Handle_Axis2.Position = [0.05, 0.08, 0.7592, 0.3];

%% Save Figures
Save_Flag = 1;
if Save_Flag == 1  
    Handle_Figure1.Color = 'none';
    set(Handle_Figure1, 'InvertHardCopy', 'off');
    print(Handle_Figure1, 'Heat_Map', '-dpdf', '-painters');
    Handle_Figure1.Color = [1, 1, 1];
    
    Handle_Figure2.Color = 'none';
    set(Handle_Figure2, 'InvertHardCopy', 'off');
    print(Handle_Figure2, 'Nominal1', '-dpdf', '-painters');
    Handle_Figure2.Color = [1, 1, 1];   
    
    Handle_Figure3.Color = 'none';
    set(Handle_Figure3, 'InvertHardCopy', 'off');
    print(Handle_Figure3, 'Nominal2', '-dpdf', '-painters');
    Handle_Figure3.Color = [1, 1, 1];
    
    Handle_Figure4.Color = 'none';
    set(Handle_Figure4, 'InvertHardCopy', 'off');
    print(Handle_Figure4, 'Nominal3', '-dpdf', '-painters');
    Handle_Figure4.Color = [1, 1, 1];
end
