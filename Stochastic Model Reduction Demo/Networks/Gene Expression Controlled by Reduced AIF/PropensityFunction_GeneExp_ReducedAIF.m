function Prop = PropensityFunction_GeneExp_ReducedAIF(x, Parameters)
% Propensity Function for Gene Expression Process Controlled by the Reduced AIF
% Controller.
% 	 Species: 		 X = [X_1; X_2; Z]
% 	 Reactions: 	R1:		 0                  --> 	X_1		        [k_0 + k*max(Z,0)]
% 				    R2:		 X_1				--> 	X_1 + X_2		[k_1*X_1]
% 				    R3:		 X_1				--> 	0				[gamma_1*X_1]
% 				    R4:		 X_2				--> 	0				[gamma_2*X_2]
%                   R5:		 0                 	-->     Z               [mu + delta_2*max(-Z,0)]
%                   R6:      Z                  -->     0               [theta*X_2 + delta_1*max(Z,0)]

%% Extract Parameters
k_0 = Parameters.k_0;
k = Parameters.k;
k_1 = Parameters.k_1;
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
mu = Parameters.mu;
theta = Parameters.theta;
delta_1 = Parameters.delta_1;
delta_2 = Parameters.delta_2;

%% Extract State Variables
X_1 = x(1);
X_2 = x(2);
Z = x(3);

%% Propensities
Prop = [ ...
        k_0 + k*max(Z,0); ...
        k_1*X_1; ...
        gamma_1*X_1; ...
		gamma_2*X_2; ...
        mu + delta_2*max(-Z,0); ...
        theta*X_2 + delta_1*max(Z,0); ...
	];
end