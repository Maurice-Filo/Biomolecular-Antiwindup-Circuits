function S = StoichiometryMatrix_PID3_AW3_Realization()
% Stoichiometry Matrix for Realized Anti Windup Topology III with 3rd order
% PID
% 	 Species:           X = [Z_1; Z_2; V_1; V_2; W_1; W_2; Z_3]
% 	 Reactions:         R1:     0                   -->         X_1                 [h_a(z_1)]
%                       R2:		0                 	-->      	Z_1                 [mu * h_1(v_1, v_2)]
%                       R3:		0                 	-->         Z_2                 [theta*h_s(y) * h_2(w_1, w_2)]
%                       R4:		Z_1 + Z_2			-->         0                   [eta * z_1 * z_2]
%                       R5:     0                   -->         V_1                 [v_0]
%                       R6:     0                   -->         V_2                 [g_1(z_1)]
%                       R7:		V_1 + V_2           -->         0                   [eta_v * v_1 * v_2]
%                       R8:     0                   -->         W_1                 [w_0]
%                       R9:     0                   -->         W_2                 [g_2(z_2)]
%                       R10:    W_1 + W_2           -->         0                   [eta_w * w_1 * w_2]
%                       R11:    Z_1                 -->         0                   [delta_1*z_1]
%                       R12:    Z_2                 -->         0                   [delta_2*z_2]
%                       R13:    V_1                 -->         0                   [delta_v1*v_1]
%                       R14:    V_2                 -->         0                   [delta_v2*v2]
%                       R15:    W_1                 -->         0                   [delta_w1*w1]
%                       R16:    W_2                 -->         0                   [delta_w2*w2]
% 				        R17:	0					-->         Z_3                 [alpha_0 / (1 + y/kappa_0)]
% 				        R18:	Z_3					-->         0                   [gamma_0*z_3]   
% 				        R19:	X_1					-->         0                   [(delta*y + delta_0*z_3) * (X_1/kappa) / (1 + X_1/kappa)]
S = [ ...
        0,      1,      0,     -1,      0,      0,      0,      0,      0,      0,     -1,      0,      0,      0,      0,      0      0,      0,      0; ...
        0,      0,      1,     -1,      0,      0,      0,      0,      0,      0,      0,     -1,      0,      0,      0,      0      0,      0,      0; ...
        0,      0,      0,      0,      1,      0,     -1,      0,      0,      0,      0,      0,     -1,      0,      0,      0      0,      0,      0; ...
        0,      0,      0,      0,      0,      1,     -1,      0,      0,      0,      0,      0,      0,     -1,      0,      0      0,      0,      0; ...
        0,      0,      0,      0,      0,      0,      0,      1,      0,     -1,      0,      0,      0,      0,     -1,      0      0,      0,      0; ...
        0,      0,      0,      0,      0,      0,      0,      0,      1,     -1,      0,      0,      0,      0,      0,     -1      0,      0,      0; ...
        0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0      1,     -1,      0; ...
	];

end

