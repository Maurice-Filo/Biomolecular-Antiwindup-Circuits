function S = StoichiometryMatrix_AIF_Reduced()
% Stoichiometry Matrix for the Reduced AIF
% 	 Species:           X = [Z]
% 	 Reactions:         R1:     0                   -->         X_1                 [h_a(max(z,0))]
%                       R2:		0                 	-->      	Z                   [mu]
%                       R3:		Z                 	-->         0                   [h_s(x_L)]
S = [ ...
        0,      1,      -1; ...
	];

end

