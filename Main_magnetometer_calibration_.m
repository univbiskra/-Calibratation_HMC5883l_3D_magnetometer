 
clc
clear
close all

% Import raw magnetometer readings from a data file ( 'my_data2.txt' my data file)

dataFileName = 'my_data2.txt'; %my_data2.txt good 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calibrationSwitch = 0; 

% Calibration Switch has two state:
% SWITCH = 0: Calibrate to unit circle (normalization).
% SWITCH = 1: Do not change the radius (scale only).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract uncalibrated magnetometer data
rawData  = importdata(dataFileName);
xUncalibrated = rawData(:,1); 
yUncalibrated = rawData(:,2); 
zUncalibrated = rawData(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create arrays to store calibrated data
XCalibrated = zeros(length(xUncalibrated),1); 
YCalibrated = zeros(length(xUncalibrated),1); 
ZCalibrated = zeros(length(xUncalibrated),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit an ellipsoid to the raw data:
% The ellipsoid equation is ax^2 + by^2 + cz^2 + 2fyz + 2gxz + 2hxy + 2px + 2qy + 2rz + d = 0.
% The vector v = [a, b, c, f, g, h, p, q, r, d]' (where k = -d).
ellipsoidParams = fit_ellipsoid(xUncalibrated,yUncalibrated,zUncalibrated); 

% Create the matrix ellipsoidMatrix (Q) and vector offsetVector (u) from the ellipsoid parameters

ellipsoidMatrix = [ellipsoidParams(1),ellipsoidParams(6),ellipsoidParams(5);ellipsoidParams(6),ellipsoidParams(2),ellipsoidParams(4);ellipsoidParams(5),ellipsoidParams(4),ellipsoidParams(3)];
offsetVector = [ellipsoidParams(7),ellipsoidParams(8),ellipsoidParams(9)]';
kValue = ellipsoidParams(10);

 % Determine calibration scaling factor and radius 
if calibrationSwitch == 1
     calibrationScale = 1;  % Do not change the radius
     calibrationRadius = sqrt(offsetVector'*(ellipsoidMatrix\offsetVector)-kValue);
elseif calibrationSwitch == 0
    calibrationScale = (1/sqrt(offsetVector'*(ellipsoidMatrix\offsetVector)-kValue)); % Normalize to unit circle
    calibrationRadius = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the calibrated offset vector offsetVectorCalibrated (b) and
% calibration matrix calibrationMatrix(Ainv)
offsetVectorCalibrated = - ellipsoidMatrix \ offsetVector; % Eqn(21)
calibrationMatrix = real(calibrationScale*sqrt(ellipsoidMatrix)); % Eqn(13) and Eqn(14) 
    
% Calibration Process:
% For each data point, the following steps are performed:
%   1. Sensor data from the raw measurements is extracted.
%   2. Calibration is applied to the data using the calibration matrix.
%   3. Calibrated values are stored in XCalibrated, YCalibrated, and ZCalibrated.

for dataIter = 1:length(xUncalibrated)
    % Sensor data
    uncalibratedData = [xUncalibrated(dataIter); yUncalibrated(dataIter); zUncalibrated(dataIter)]; 
    
    % Apply calibration matrix and offset vector
    calibratedData = calibrationMatrix*(uncalibratedData - offsetVectorCalibrated);
    
    % Store calibrated values
    XCalibrated(dataIter) = calibratedData(1);
    YCalibrated(dataIter) = calibratedData(2);
    ZCalibrated(dataIter) = calibratedData(3);
end
    
% Plot uncalibrated data
subplot(1,2,1);
plot_fitted_ellipsoid(ellipsoidParams); 
hold on;

scatter3(xUncalibrated,yUncalibrated,zUncalibrated,'fill','MarkerFaceColor','red');
title({'Before magnetometer calibration','(Ellipsoid fitted)'});
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
axis equal;

% Plot calibrated data
subplot(1,2,2);
plot_sphere_shape([0,0,0],calibrationRadius);
hold on;

scatter3(XCalibrated,YCalibrated,ZCalibrated,'fill','MarkerFaceColor','blue');
if calibrationSwitch == 0
    title({'After magnetometer calibration','(Normalized to unit circle)'});
else
    title({'After magnetometer calibration'});
end
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
axis equal;

% Print calibration matrices
fprintf('3D magnetometer calibration via ellipsoid fitting');
fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
fprintf('\nThe Calibration Equation to be Applied:') 
fprintf('\n\t\t\t\th_calibrated = CalibrationMatrix * (h_uncalibrated - OffsetVector) \nWhere,')
fprintf('\nh_uncalibrated   = Measured sensor data vector');
fprintf('\nh_hat = Vector of Calibrated Sensor Data');
fprintf('\nh_calibrated = Vector of Calibrated Sensor Data');
fprintf('\n\nCalibrationMatrix =\n'); disp(calibrationMatrix);
fprintf('\nOffsetVector =\n'); disp(offsetVectorCalibrated);
 
 