function Propensity = PropensityFunction_AW3_Reduced(z, y, Parameters)
% Propensity Function for the Reduced Anti-Windup Topology III
% 	 Species:           X = [Z]
% 	 Reactions:         R1:     0                   -->         X_1                 [k*h_a(max(z,0))]
%                       R2:		0                 	-->      	Z                   [mu * h_1(max(z,0))]
%                       R3:		Z                 	-->         0                   [theta*h_s(x_L)*h_2(max(-z,0))]

%% Extract Parameters
mu = Parameters.mu;

%% Construct Propensity Function
Propensity = [  Parameters.ActuationFunction(max(z, 0), Parameters.ActuationParameters); ...
                mu * Parameters.h_1(max(z, 0), Parameters.WindupParameters); ...
                Parameters.SensingFunction(y, Parameters.SensingParameters) * Parameters.h_2(max(-z,0), Parameters.WindupParameters); ...
              ];
end

