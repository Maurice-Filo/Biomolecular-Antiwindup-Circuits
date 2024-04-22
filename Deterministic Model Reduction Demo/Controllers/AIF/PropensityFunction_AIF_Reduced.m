function Propensity = PropensityFunction_AIF_Reduced(z, y, Parameters)
% Propensity Function for the Reduced AIF
% 	 Species:           X = [Z]
% 	 Reactions:         R1:     0                   -->         X_1                 [h_a(max(z,0))]
%                       R2:		0                 	-->      	Z                   [mu]
%                       R3:		Z                 	-->         0                   [h_s(x_L)]

%% Extract Parameters
mu = Parameters.mu;

%% Construct Propensity Function
Propensity = [  Parameters.ActuationFunction(max(z, 0), Parameters.ActuationParameters); ...
                mu; ...
                Parameters.SensingFunction(y, Parameters.SensingParameters); ...
              ];
end

