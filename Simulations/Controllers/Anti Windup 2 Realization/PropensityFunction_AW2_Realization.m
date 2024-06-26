function Propensity = PropensityFunction_AW2_Realization(z, y, Parameters)
% Propensity Function for Realized Anti Windup Topology II
% 	 Species:           X = [Z_1; Z_2; V_1; V_2; W_1; W_2]
% 	 Reactions:         R1:     0                   -->         X_1                 [h_a(z_1)]
%                       R2:		0                 	-->      	Z_1                 [mu]
%                       R3:		0                 	-->         Z_2                 [theta*h_s(y)]
%                       R4:		Z_1 + Z_2			-->         0                   [eta * z_1 * z_2]
%                       R5:     0                   -->         V_1                 [v_0]
%                       R6:     0                   -->         V_2                 [g_1(z_1)]
%                       R7:		V_1 + V_2           -->         0                   [eta_v * v_1 * v_2]
%                       R8:     0                   -->         W_1                 [w_0]
%                       R9:     0                   -->         W_2                 [g_2(z_2)]
%                       R10:    W_1 + W_2           -->         0                   [eta_w * w_1 * w_2]
%                       R11:    Z_1                 -->         0                   [(delta_1 + h_1(v_1, v_2)) * z_1]
%                       R12:    Z_2                 -->         0                   [(delta_2 + h_2(w_1, w_2)) * z_2]
%                       R13:    V_1                 -->         0                   [delta_v1*v_1]
%                       R14:    V_2                 -->         0                   [delta_v2*v2]
%                       R15:    W_1                 -->         0                   [delta_w1*w1]
%                       R16:    W_2                 -->         0                   [delta_w2*w2]

%% Extract Parameters
mu = Parameters.mu;
eta = Parameters.eta;
eta_v = Parameters.eta_v;
eta_w = Parameters.eta_w;
v_0 = Parameters.v_0;
w_0 = Parameters.w_0;
delta_1 = Parameters.delta_1;
delta_2 = Parameters.delta_2;
delta_v1 = Parameters.delta_v1;
delta_v2 = Parameters.delta_v2;
delta_w1 = Parameters.delta_w1;
delta_w2 = Parameters.delta_w2;

%% Extract Functions
h_a = Parameters.ActuationFunction;
ActuationParameters = Parameters.ActuationParameters;
h_s = Parameters.SensingFunction;
SensingParameters = Parameters.SensingParameters;
h_1 = Parameters.h_1;
h_2 = Parameters.h_2;
g_1 = Parameters.g_1;
g_2 = Parameters.g_2;

%% Extract State Variables
z_1 = z(1);
z_2 = z(2);
v_1 = z(3);
v_2 = z(4);
w_1 = z(5);
w_2 = z(6);

%% Construct Propensity Function
Propensity = [  h_a(z_1, ActuationParameters); ...
                mu; ...
                h_s(y, SensingParameters); ...
                eta * z_1 * z_2; ...
                v_0; ...
                g_1(z_1); ...
                eta_v * v_1 * v_2; ...
              	w_0; ...
                g_2(z_2); ...
                eta_w * w_1 * w_2; ...
                (delta_1 + h_1(v_1, v_2)) * z_1; ...
                (delta_2 + h_2(w_1, w_2)) * z_2; ...
                delta_v1 * v_1; ...
                delta_v2 * v_2; ...
                delta_w1 * w_1; ...
                delta_w2 * w_2; ...
              ];
end

