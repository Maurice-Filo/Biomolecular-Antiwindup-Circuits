function y = PlasmaInsulin(time) %#codegen
%PlasmaInsulin - Helper function for SimBiology examples.
%
% y = PlasmaInsulin(time);
% where time is non-dimensionalized by scaling by 1/hour, and y is
% non-dimensionalized by scaling by liter/picomole.

%   Copyright 2011-2021 The MathWorks, Inc.

% time non-dimensional, and scaled by 1/hour
% y is the non-dimensionalized plasma insulin concentration. It's multiplied by the unit
% picomole/liter.

% timeData and yData were obtained from Figure 1 of Meal Simulation Model
% of the Glucose-Insulin System. C. Dalla Man, R.A. Rizza, and C. Cobelli.
% IEEE Transactions on Biomedical Engineering (2007) 54(10), 1740-1749.

% time data in hour
timeData = [0; 0.09; 0.15; 0.23; 0.33; 0.49; 0.66; 0.82; 0.99; 1.25; ...
    1.50; 1.98; 2.49; 2.98; 3.48; 3.99; 4.31; 4.64; 4.99; 5.98; 7];

% plamsa insulin data in picomole/liter
yData = [24; 30; 45; 112; 184; 343; 308; 317; 355; 309; 335; 252; ...
    180; 116; 83; 76; 62; 56; 57; 41; 31];

% Use linear interpolation
y = interp1(timeData, yData, real(time));
