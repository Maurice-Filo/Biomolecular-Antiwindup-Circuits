function Parameters = Parameters_GeneExp_AIF()
% Initial Parameters for Gene Expression Process Controlled by the AIF
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

Parameters.k_0 = 0;
Parameters.k = 1;
Parameters.k_0 = 10;
Parameters.k_1 = 1;
Parameters.gamma_1 = 5;
Parameters.gamma_2 = 5;
Parameters.mu = 10;
Parameters.theta = 1;
Parameters.eta = 1e4;
Parameters.delta_1 = 0;
Parameters.delta_2 = 0;
end

