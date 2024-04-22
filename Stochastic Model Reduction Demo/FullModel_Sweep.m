%% Clear Workspace
close all
clear
clc
Cluster = 0;

%% Configure Working Environment
if Cluster == 1
    ProjectPath = '/cluster/home/mfilo/Code';
    DataPath = fullfile(ProjectPath, 'Data');
    addpath(genpath(ProjectPath));
    euler = parcluster('local');
    pool = parpool(euler, 24);
else
    ProjectPath = '/Users/mfilo/Library/CloudStorage/GoogleDrive-maurice.g.filo@gmail.com/My Drive/Research at ETH/My Papers and Presentations/Protection Circuits/Matlab/Stochastic Model Reduction Demo';
    DataPath = fullfile(ProjectPath, 'Data');
    addpath(genpath(ProjectPath));
end

%% Network Construction
StoichiometryMatrix = StoichiometryMatrix_GeneExp_AIF();
PropensityFunction = @PropensityFunction_GeneExp_AIF;
Parameters = Parameters_GeneExp_AIF();

%% Simulation Settings
IC = [0; 0; 0; 0];
TimeSpan = [0, 40];
N_Trajectories = 1e3;
N_TimeGrid = 100;
N_ParameterGrid = 2;

%% Swept Parameters
gamma_1_vector = linspace(0.5, 5, N_ParameterGrid);
gamma_2_vector = linspace(0.5, 5, N_ParameterGrid);

%% Moments Computations
SimTime = tic;
Mean = cell(N_ParameterGrid, N_ParameterGrid);
Covariance = cell(N_ParameterGrid, N_ParameterGrid);
Skewness = cell(N_ParameterGrid, N_ParameterGrid);
for i = 1 : N_ParameterGrid
    Parameters.gamma_1 = gamma_1_vector(i);
    parfor j = 1 : N_ParameterGrid
        [i,j]
        localParameters = Parameters; 
        localParameters.gamma_2 = gamma_2_vector(j);
        [~, Mean{i,j}, Covariance{i,j}, ~, Skewness{i,j}, ~] = GenerateMoments(StoichiometryMatrix, PropensityFunction, localParameters, IC, TimeSpan, N_Trajectories, N_TimeGrid);
    end
end
T = linspace(TimeSpan(1), TimeSpan(2), N_TimeGrid);
SimTime = toc(SimTime);

%% Save data
Save_Flag = 1;
FileName = fullfile(DataPath, 'FullModel');
if Save_Flag == 1
    clear euler pool
    save(FileName);
end