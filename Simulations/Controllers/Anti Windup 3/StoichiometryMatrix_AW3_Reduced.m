function S = StoichiometryMatrix_AW3_Reduced()
% Stoichiometry Matrix for the Reduced Anti-Windup Topology III
% 	 Species:           X = [Z]
% 	 Reactions:         R1:     0                   -->         X_1                 [k*h_a(max(z,0))]
%                       R2:		0                 	-->      	Z                   [mu * h_1(max(z,0))]
%                       R3:		Z                 	-->         0                   [theta*h_s(x_L)*h_2(max(-z,0))]
S = [ ...
        0,      1,      -1; ...
	];

end

