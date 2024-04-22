function Propensity = PropensityFunction_GeneExp(x, Parameters)
% Propensity Function for Gene Expression Network
% 	 Species:           X = [X_1; X_2]
% 	 Reactions:         R1:      0                      -->         X_1                     [k_0]  
%                       R2:		 X_1					-->         X_1 +  X_2				[k_1*x_1]
%                       R3:		 X_1					-->         0                       [gamma_1*x_1]
%                       R4:		 X_2					-->         0                       [gamma_2*x_2]

%% Extract Parameters
k_0 = Parameters.k_0;
k_1 = Parameters.k_1;
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;

%% Construct Propensity Function
Propensity = [  k_0; ...
                k_1 * x(1); ...
                gamma_1 * x(1); ...
                gamma_2 * x(2); ...
              ];
end

