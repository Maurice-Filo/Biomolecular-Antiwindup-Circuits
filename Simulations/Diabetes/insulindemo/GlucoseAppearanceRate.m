function y = GlucoseAppearanceRate(time) %#codegen
%GlucoseAppearanceRate - Helper function for SimBiology examples.
% y = GlucoseAppearanceRate(time);

%   Copyright 2011-2021 The MathWorks, Inc.

% time is non-dimensionalized by scaling by 1/hour
% y is the non-dimensionalized glucose rate of appearance scaled by minute/milligram

% timeData and yData were obtained from Figure 1 of Meal Simulation Model
% of the Glucose-Insulin System. C. Dalla Man, R.A. Rizza, and C. Cobelli.
% IEEE Transactions on Biomedical Engineering (2007) 54(10), 1740-1749.

% time data in hour
timeData = [0; 0.21; 0.29; 0.42; 0.57; 0.73; 0.91; 1.13; 1.38; 1.74; ...
    2.25; 2.76; 3.25; 3.74; 4.17; 4.50; 4.83; 5.50; 6.51; 7.00];

% glucose appearance in milligram/minute
yData = [0; 9.13; 9.22; 9.75; 7.49; 7.21; 6.77; 6.27; 4.87; ...
    4.43; 3.31; 2.54; 2.29; 1.84; 1.83; 1.31; 1.26; 1.12; 0.54; 0.54];

% Use linear interpolation
y = interp1(timeData, yData, real(time));
