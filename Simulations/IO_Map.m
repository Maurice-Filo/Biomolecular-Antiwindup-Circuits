%% Clear Workspace
close all
clear
clc

Factor = 1e3;
% Plant Selection
StoichiometryMatrix_Plant = StoichiometryMatrix_Star();
PropensityFunction_Plant = @PropensityFunction_Star;
% Parameters
Parameters_Plant.b_0 = .2;
Parameters_Plant.b_1 = .2;
Parameters_Plant.b_2 = .2;
Parameters_Plant.b_3 = .2;
Parameters_Plant.b_4 = .2;
Parameters_Plant.b_5 = .2;
Parameters_Plant.k_1 = 0.1 * Factor;
Parameters_Plant.k_2 = 0.1 * Factor;
Parameters_Plant.k_3 = 0.1 * Factor;
Parameters_Plant.k_4 = 0.1 * Factor;
Parameters_Plant.k_5 = 0.1 * Factor;
Parameters_Plant.kappa_1 = Factor;
Parameters_Plant.kappa_2 = Factor;
Parameters_Plant.kappa_3 = Factor;
Parameters_Plant.kappa_4 = Factor;
Parameters_Plant.kappa_5 = Factor;
Parameters_Plant.gamma_F = 0.3;
Parameters_Plant.kappa_F = 1e-2;
Parameters_Plant.gamma_1 = 0.1;
Parameters_Plant.gamma_2 = 0.1;
Parameters_Plant.gamma_3 = 0.1;
Parameters_Plant.gamma_4 = 0.1;
Parameters_Plant.gamma_5 = 0.1;
Parameters_Plant.gamma_6 = 0.1;
L = 6;                  % Number of Plant Species
Input_Index = 1;
Output_Index = 6;

%% Input Range
u = 0 : 0.1 : 0.06 * Factor;

%% Output 
y = 0 * u;
for i = 1 : length(u)
    Parameters_Plant.b_0 = u(i);
    Fun = @(x) RHS(x, StoichiometryMatrix_Plant, PropensityFunction_Plant, Parameters_Plant);
    X = fsolve(Fun, zeros(L,1));
    y(i) = X(Output_Index);
end

%% Plot
plot(u, y);
    
function dx = RHS(x, StoichiometryMatrix, PropensityFunction, Parameters_Plant)
    dx = StoichiometryMatrix * PropensityFunction(x, Parameters_Plant);
end