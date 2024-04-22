function Y = ResampleZOH(X, T, t_vector)
% ResampleZOH performs zero-order hold (ZOH) resampling of a sequence.
%
% Inputs:
% - X: Data matrix to be resampled. Each column is a data point and rows are dimensions/variables.
% - T: Time vector associated with the data points in X.
% - t_vector: New time vector for which resampled values are needed.
%
% Outputs:
% - Y: Resampled data using ZOH.

%% Initialize output matrix
Y = zeros(size(X,1),length(t_vector));
    j_Previous = 1;
    for i = 1 : length(t_vector)
        t_i = t_vector(i);
        for j = j_Previous : length(T)-1
            if (t_i >= T(j)) && (t_i < T(j+1))
                Y(:,i) = X(:,j);
                j_Previous = j;
                break;
            end
        end
    end
end

