function StochasticPlot_Statistics(T, X, t_vector, MeanX, VarX, Parameters, Title, Names, N_DisplayedTrajectories)
    %% Create Colors for Plant and Controller Species
    Colors = lines(2);
    %% Total Number of Species
    N = size(X{1},1);
    %% Create Figure and Subplots 
    Handle_Figure = figure();
    Handle_SamplePath = gobjects(N,1);
    Handle_Mean = gobjects(N,1);
    Handle_Variance = gobjects(N,1);
    Handle_Distribution = gobjects(N,1);
    for k = 1 : N
        Handle_SamplePath(k) = subplot(N, 4, 4*(k-1) + 1);
        Handle_Mean(k) = subplot(N, 4, 4*(k-1) + 2);
        Handle_Variance(k) = subplot(N, 4, 4*(k-1) + 3);
        Handle_Distribution(k) = subplot(N, 4, 4*(k-1) + 4);
    end
    %% Main Title
    sgtitle([Title, ', Stochastic Setting']);
    %% Plotting...
    for k = 1 : N
        if N-k < 2
            c = Colors(1,:);
        else 
            c = Colors(2,:);
        end
        % Sample Paths
        hold(Handle_SamplePath(k), 'on');
        for i = 1 : N_DisplayedTrajectories
            stairs(Handle_SamplePath(k), T{i}, X{i}(k,:), 'Color', [c, 0.2], 'LineWidth', Parameters.LineWidth);
        end
        Handle_SamplePath(k).YLabel.String = ['$', Names{k}, '$'];
        Handle_SamplePath(k).YLabel.Interpreter = 'latex';
        Handle_SamplePath(k).XGrid = 'on';
        Handle_SamplePath(k).YGrid = 'on';
        Handle_SamplePath(k).XMinorGrid = 'on';
        Handle_SamplePath(k).YMinorGrid = 'on';
        Handle_SamplePath(k).FontSize = Parameters.FontSize;
        Handle_SamplePath(k).FontName = 'Times New Roman';
        
        % Means
        hold(Handle_Mean(k), 'on');
        plot(Handle_Mean(k), t_vector, MeanX(k,:), 'Color', c);
        Handle_Mean(k).XGrid = 'on';
        Handle_Mean(k).YGrid = 'on';
        Handle_Mean(k).XMinorGrid = 'on';
        Handle_Mean(k).YMinorGrid = 'on';
        Handle_Mean(k).FontSize = Parameters.FontSize;
        Handle_Mean(k).FontName = 'Times New Roman';
        
        % Variances
        hold(Handle_Variance(k), 'on');
        plot(Handle_Variance(k), t_vector, VarX(k,:), 'Color', c);
        Handle_Variance(k).XGrid = 'on';
        Handle_Variance(k).YGrid = 'on';
        Handle_Variance(k).XMinorGrid = 'on';
        Handle_Variance(k).YMinorGrid = 'on';
        Handle_Variance(k).FontSize = Parameters.FontSize;
        Handle_Variance(k).FontName = 'Times New Roman';
        
        % Final Time Distribution
        Xk_bar = zeros(length(T),1);
        for i = 1 : length(T)
            Xk_bar(i) = X{i}(k,end);
        end
        hold(Handle_Distribution(k), 'on');
        histogram(Handle_Distribution(k), Xk_bar, 'FaceColor', c);
        Handle_Distribution(k).XGrid = 'on';
        Handle_Distribution(k).YGrid = 'on';
        Handle_Distribution(k).XMinorGrid = 'on';
    	Handle_Distribution(k).YMinorGrid = 'on';
    	Handle_Distribution(k).FontSize = Parameters.FontSize;
    	Handle_Distribution(k).FontName = 'Times New Roman';
    end 
    Handle_SamplePath(1).Title.String = 'Trajectories';
    Handle_Mean(1).Title.String = 'Expectation';
    Handle_Variance(1).Title.String = 'Variance';
    Handle_Distribution(1).Title.String = 'Final Distribution';
    Handle_SamplePath(N).XLabel.String = '$t$';
    Handle_SamplePath(N).XLabel.Interpreter = 'latex';
    Handle_Mean(N).XLabel.String = '$t$';
    Handle_Mean(N).XLabel.Interpreter = 'latex';
    Handle_Variance(N).XLabel.String = '$t$';
    Handle_Variance(N).XLabel.Interpreter = 'latex';
end

