function S = StoichiometryMatrix_GeneExp()
% Stoichiometry Matrix for Gene Expression Network
% 	 Species:           X = [X_1; X_2]
% 	 Reactions:         R1:      0                      -->         X_1                     [k_0]  
%                       R2:		 X_1					-->         X_1 +  X_2				[k_1*x_1]
%                       R3:		 X_1					-->         0                       [gamma_1*x_1]
%                       R4:		 X_2					-->         0                       [gamma_2*x_2]
S = [ ...
        1,      0,	   -1,		0; ...
		0,      1,		0,	   -1; ...
	];

end

