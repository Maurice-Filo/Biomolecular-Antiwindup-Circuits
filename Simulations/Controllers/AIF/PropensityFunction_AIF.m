function Propensity = PropensityFunction_AIF(z, y, Parameters)
% Propensity Function for AIF
% 	 Species:           X = [Z_1; Z_2]
% 	 Reactions:         R1:     0                   -->         X_1                 [h_a(z_1)]
%                       R2:		0                 	-->      	Z_1                 [mu]
%                       R3:		0                 	-->         Z_2                 [h_s(x_L)]
%                       R4:		Z_1 + Z_2			-->         0                   [eta * z_1 * z_2]
%                       R5:     Z_1                 -->         0                   [delta_1*z_1]
%                       R6:     Z_2                 -->         0                   [delta_2*z_2]

%% Extract Parameters
mu = Parameters.mu;
eta = Parameters.eta;
delta_1 = Parameters.delta_1;
delta_2 = Parameters.delta_2;

%% Construct Propensity Function
Propensity = [  Parameters.ActuationFunction(z(1), Parameters.ActuationParameters); ...
                mu; ...
                Parameters.SensingFunction(y, Parameters.SensingParameters); ...
                eta * z(1) * z(2); ...
                delta_1 * z(1); ...
                delta_2 * z(2); ...
              ];
end

