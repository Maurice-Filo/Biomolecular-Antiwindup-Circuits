function S = StoichiometryMatrix_Sequestration()
% Stoichiometry Matrix for Sequestration
% 	 Species:           X = [V_1; V_2]
% 	 Reactions:         R1:		0                 	-->      	V_1                 [u_0]
%                       R2:		0                 	-->         V_2                 [u]
%                       R3:		V_1 + V_2			-->         0                   [eta * v_1 * v_2]
%                       R4:     V_1                 -->         0                   [delta_1*v_1]
%                       R5:     V_2                 -->         0                   [delta_2*v_2]
S = [ ...
        1,      0,     -1,     -1,      0; ...
        0,      1,     -1,      0,     -1; ...
	];

end

