function Propensity = PropensityFunction_Sequestration(v, Parameters, t)
% Propensity Function for Sequestration
% 	 Species:           X = [V_1; V_2]
% 	 Reactions:         R1:		0                 	-->      	V_1                 [u_0]
%                       R2:		0                 	-->         V_2                 [u]
%                       R3:		V_1 + V_2			-->         0                   [eta * v_1 * v_2]
%                       R4:     V_1                 -->         0                   [delta_1*v_1]
%                       R5:     V_2                 -->         0                   [delta_2*v_2]

%% Extract Parameters
u = Parameters.u;
eta = Parameters.eta;
u_0 = Parameters.u_0;
delta_1 = Parameters.delta_1;
delta_2 = Parameters.delta_2;

%% Construct Propensity Function
Propensity = [  u_0; ...
                u(t); ...
                eta * v(1) * v(2); ...
                delta_1 * v(1); ...
                delta_2 * v(2); ...
              ];
end

