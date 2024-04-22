function a = PropensityFunction_Star(x, Parameters)
% Propensity Function for Star Network
% 	 Species: 		X = [X_1; X_2; X_3; X_4; X_5; X_6]
% 	 Reactions: 	
%                   R1:     0           -->     X_1                 [b_0]
% 				    R2:     X_1			--> 	X_1 +  X_2          [k_1*X_1/(X_1 + kappa_1) + b_1]
% 				    R3:		X_2			--> 	X_2 +  X_3          [k_2*X_2/(X_2 + kappa_2) + b_2]
% 				    R4:		X_3			--> 	X_3 +  X_4          [k_3*X_3/(X_3 + kappa_3) + b_3]
% 				    R5:     X_4			--> 	X_4 +  X_5          [k_4*X_4/(X_4 + kappa_4) + b_4]
% 				    R6:     X_5			--> 	X_5 +  X_6          [k_5*X_5/(X_5 + kappa_5) + b_5]
% 				    R7:     X_2			--> 	0                   [gamma_F*X_6*X_2 / (X_2 + kappa_F)]
%                   R8:		X_1			--> 	0                   [gamma_1*X_1]
% 				    R9:		X_2			--> 	0                   [gamma_2*X_2]
% 				    R10:	X_3			--> 	0                   [gamma_3*X_3]
% 				    R11:    X_4			--> 	0                   [gamma_4*X_4]
% 				    R12:	X_5			--> 	0                   [gamma_5*X_5]
% 				    R13:	X_6			--> 	0                   [gamma_6*X_6]

%% Extract Parameters
gamma_F = Parameters.gamma_F;
kappa_F = Parameters.kappa_F;
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
gamma_3 = Parameters.gamma_3;
gamma_4 = Parameters.gamma_4;
gamma_5 = Parameters.gamma_5;
gamma_6 = Parameters.gamma_6;
k_1 = Parameters.k_1;
k_2 = Parameters.k_2;
k_3 = Parameters.k_3;
k_4 = Parameters.k_4;
k_5 = Parameters.k_5;
b_0 = Parameters.b_0;
b_1 = Parameters.b_1;
b_2 = Parameters.b_2;
b_3 = Parameters.b_3;
b_4 = Parameters.b_4;
b_5 = Parameters.b_5;
kappa_1 = Parameters.kappa_1;
kappa_2 = Parameters.kappa_2;
kappa_3 = Parameters.kappa_3;
kappa_4 = Parameters.kappa_4;
kappa_5 = Parameters.kappa_5;

%% Extract State Variables
X_1 = x(1);
X_2 = x(2);
X_3 = x(3);
X_4 = x(4);
X_5 = x(5);
X_6 = x(6);

%% Propensities
a = [ ...
        b_0; ...
        b_1 + k_1*X_1/(X_1 + kappa_1); ... 
        b_2 + k_2*X_2/(X_2 + kappa_2); ... 
        b_3 + k_3*X_3/(X_3 + kappa_3); ... 
        b_4 + k_4*X_4/(X_4 + kappa_4); ... 
        b_5 + k_5*X_5/(X_5 + kappa_5); ... 
		gamma_F*X_6*X_2 / (X_2 + kappa_F); ...
		gamma_1*X_1; ...
		gamma_2*X_2; ...
		gamma_3*X_3; ...
		gamma_4*X_4; ...
		gamma_5*X_5; ...
		gamma_6*X_6; ...
	];
end