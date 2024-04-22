%% Clear Workspace
close all
clear
clc

%% Load Results
Data_Full = load('FullModel.mat');
Data_Reduced = load('ReducedModel.mat');

%% Nominal Simulations
Nominal_Indeces = [2, 2];%; ...
                   %1, 2; ...
                   %2, 1; ...
                   %2, 2];

%% Plot Results
Colors = lines(10);
PlotParameters.FontSize = 24;
PlotParameters.LineWidth = 3;
Title = 'Gene Expression Process';
Names = {'X_1', 'X_2', 'Z_1', 'Z_2'};
N_Species = size(Data_Full.StoichiometryMatrix, 1);
Color_Full = Colors(1,:);
Color_Reduced = Colors(2, :);

for counter = 1 : size(Nominal_Indeces, 1)
    i = Nominal_Indeces(counter, 1);
    j = Nominal_Indeces(counter, 2);
    % Create Figure and Subplots for the Mean 
    figure('Position', [100, 100, 1000, 500]);
    Handle_Mean = gobjects(N_Species,1);
    for k = 1 : N_Species
	    Handle_Mean(k) = subplot(N_Species, 1, k);
        Handle_Mean(k).XGrid = 'on';
	    Handle_Mean(k).YGrid = 'on';
	    Handle_Mean(k).XMinorGrid = 'on';
	    Handle_Mean(k).YMinorGrid = 'on';
	    Handle_Mean(k).FontSize = PlotParameters.FontSize;
        Handle_Mean(k).XLabel.String = 't';
        Handle_Mean(k).YLabel.String = ['$E[', Names{k}, '(t)]$'];
        Handle_Mean(k).XLabel.Interpreter = 'latex';
        Handle_Mean(k).YLabel.Interpreter = 'latex';
        sgtitle([Title, ', Expectations']);
        hold(Handle_Mean(k), 'on');
	    plot(Handle_Mean(k), Data_Full.T, Data_Full.Mean{i,j}(k,:), 'Color', Color_Full);
        if k <= 3
        plot(Handle_Mean(k), Data_Reduced.T, Data_Reduced.Mean{i,j}(k,:), 'Color', Color_Reduced);
        end
    end

    % Create Figure and Subplots for the Covariance 
    figure('Position', [100, 100, 1000, 500]);
    Handle_Covariance = gobjects(N_Species^2,1);
    for k = 1 : N_Species^2
	    Handle_Covariance(k) = subplot(N_Species, N_Species, k);
        Handle_Covariance(k).XGrid = 'on';
	    Handle_Covariance(k).YGrid = 'on';
	    Handle_Covariance(k).XMinorGrid = 'on';
	    Handle_Covariance(k).YMinorGrid = 'on';
	    Handle_Covariance(k).FontSize = PlotParameters.FontSize;
        Handle_Covariance(k).XLabel.String = 't';
        Handle_Covariance(k).YLabel.String = ['$C[', Names{floor((k - 1) / N_Species) + 1}, '(t),' Names{mod(k - 1, N_Species) + 1}, '(t)]$'];
        Handle_Covariance(k).XLabel.Interpreter = 'latex';
        Handle_Covariance(k).YLabel.Interpreter = 'latex';
        sgtitle([Title, ', Covariances']);
        hold(Handle_Covariance(k), 'on');
	    plot(Handle_Covariance(k), Data_Full.T, squeeze(Data_Full.Covariance{i,j}(floor((k - 1) / N_Species) + 1, mod(k - 1, N_Species) + 1, :)), 'Color', Color_Full);
        if (floor((k - 1) / N_Species) + 1) <= 3 &&  (mod(k - 1, N_Species) + 1 <=3)
        plot(Handle_Covariance(k), Data_Reduced.T, squeeze(Data_Reduced.Covariance{i,j}(floor((k - 1) / N_Species) + 1, mod(k - 1, N_Species) + 1, :)), 'Color', Color_Reduced);
        end
    end
end