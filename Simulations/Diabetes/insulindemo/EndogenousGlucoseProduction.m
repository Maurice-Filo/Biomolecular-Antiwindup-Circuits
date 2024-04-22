function y = EndogenousGlucoseProduction(time) %#codegen
%EndogenousGlucoseProduction - Helper function for SimBiology examples.
% y = EndogenousGlucoseProduction(time);

%   Copyright 2011-2021 The MathWorks, Inc.

% time is non-dimensionalized by scaling by 1/hour.
% y is the non-dimensionalized endogenous glucose production, scaled by minute/milligram.

% timeData and yData were obtained from Figure 1 of Meal Simulation Model
% of the Glucose-Insulin System. C. Dalla Man, R.A. Rizza, and C. Cobelli.
% IEEE Transactions on Biomedical Engineering (2007) 54(10), 1740-1749.

% time data in hour
timeData = [0; 0.04; 0.13; 0.22; 0.28; 0.40; 0.58; 0.73; 0.91; 1.11; ...
    1.36; 1.72; 2.21; 2.71; 3.21; 3.71; 4.12; 4.45; 4.76; 5.44; 6.43; 7];

% endogenous glucose production in milligram/minute

yData = [1.91; 1.67; 1.38; 1.17; 1.06; 0.73; 0.51; 0.45; 0.47; 0.51; ...
    0.53; 0.55; 0.63; 0.71; 0.74; 0.82; 0.70; 0.99; 0.98; 1.08; 1.10; 1.09];

% Use linear interpolation
y = interp1(timeData, yData, real(time));
