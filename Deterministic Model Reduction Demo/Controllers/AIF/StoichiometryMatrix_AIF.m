function S = StoichiometryMatrix_AIF()
% Stoichiometry Matrix for AIF
% 	 Species:           X = [Z_1; Z_2]
% 	 Reactions:         R1:     0                   -->         X_1                 [h_a(z_1)]
%                       R2:		0                 	-->      	Z_1                 [mu]
%                       R3:		0                 	-->         Z_2                 [h_s(x_L)]
%                       R4:		Z_1 + Z_2			-->         0                   [eta * z_1 * z_2]
%                       R5:     Z_1                 -->         0                   [delta_1*z_1]
%                       R6:     Z_2                 -->         0                   [delta_2*z_2]
S = [ ...
        0,      1,      0,     -1,     -1,      0; ...
        0,      0,      1,     -1,      0,     -1; ...
	];

end

