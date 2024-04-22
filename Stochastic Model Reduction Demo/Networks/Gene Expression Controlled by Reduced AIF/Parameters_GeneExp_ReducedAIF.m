function Parameters = Parameters_GeneExp_ReducedAIF()
% Initial Parameters for Gene Expression Process Controlled by the Reduced AIF
% Controller.
% 	 Species: 		 X = [X_1; X_2; Z]
% 	 Reactions: 	R1:		 0                  --> 	X_1		        [k_0 + k*max(Z,0)]
% 				    R2:		 X_1				--> 	X_1 + X_2		[k_1*X_1]
% 				    R3:		 X_1				--> 	0				[gamma_1*X_1]
% 				    R4:		 X_2				--> 	0				[gamma_2*X_2]
%                   R5:		 0                 	-->     Z               [mu + delta_2*max(-Z,0)]
%                   R6:      Z                  -->     0               [theta*X_2 + delta_1*max(Z,0)]

Parameters.k_0 = 0;
Parameters.k = 1;
Parameters.k_0 = 10;
Parameters.k_1 = 1;
Parameters.gamma_1 = 5;
Parameters.gamma_2 = 5;
Parameters.mu = 10;
Parameters.theta = 1;
Parameters.delta_1 = 0;
Parameters.delta_2 = 0;
end

