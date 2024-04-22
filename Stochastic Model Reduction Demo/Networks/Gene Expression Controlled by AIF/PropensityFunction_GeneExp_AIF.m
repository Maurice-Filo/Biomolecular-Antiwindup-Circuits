function Prop = PropensityFunction_GeneExp_AIF(x, Parameters)
% Propensity Function for Gene Expression Process Controlled by the AIF
% Controller.
% 	 Species: 		 X = [X_1; X_2; Z_1; Z_2]
% 	 Reactions: 	R1:		 Z_1                --> 	Z_1 + X_1		[k_0 + k*Z_1]
% 				    R2:		 X_1				--> 	X_1 + X_2		[k_1*X_1]
% 				    R3:		 X_1				--> 	0				[gamma_1*X_1]
% 				    R4:		 X_2				--> 	0				[gamma_2*X_2]
%                   R5:		 0                 	-->     Z_1             [mu]
%                   R6:		 0                 	-->     Z_2             [theta*X_2]
%                   R7:		 Z_1 + Z_2			-->     0               [eta*Z_1*Z_2]
%                   R8:      Z_1                -->     0               [delta_1*Z_1]
%                   R9:      Z_2                -->     0               [delta_2*Z_2]

%% Extract Parameters
k_0 = Parameters.k_0;
k = Parameters.k;
k_1 = Parameters.k_1;
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
mu = Parameters.mu;
theta = Parameters.theta;
eta = Parameters.eta;
delta_1 = Parameters.delta_1;
delta_2 = Parameters.delta_2;

%% Extract State Variables
X_1 = x(1);
X_2 = x(2);
Z_1 = x(3);
Z_2 = x(4);

%% Propensities
Prop = [ ...
        k_0 + k*Z_1; ...
        k_1*X_1; ...
        gamma_1*X_1; ...
		gamma_2*X_2; ...
        mu; ...
        theta*X_2; ...
        eta*Z_1*Z_2; ...
        delta_1*Z_1; ...
        delta_2*Z_2; ...
	];
end