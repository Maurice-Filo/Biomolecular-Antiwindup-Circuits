function S = StoichiometryMatrix_Star()
% Stoichiometry Matrix for Star Network
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


S = [  1,	  0,	  0,	  0,	  0,	  0,	  0,	 -1,	  0,	  0,	  0,	  0,	  0; ...
	   0,	  1,	  0,	  0,	  0,	  0,	 -1,      0,	 -1,	  0,	  0,	  0,	  0; ...
	   0,     0,	  1,	  0,	  0,	  0,	  0,      0,	  0,	 -1,	  0,	  0,	  0; ...
	   0,	  0,	  0,	  1,	  0,	  0,	  0,      0,	  0,	  0,	 -1,	  0,	  0; ...
	   0,	  0,	  0,	  0,	  1,	  0,	  0,      0,	  0,	  0,	  0,	 -1,	  0; ...
	   0,     0,	  0,	  0,	  0,	  1,	  0,      0,	  0,	  0,	  0,	  0,	 -1; ...
       ];
end