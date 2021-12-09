% Plot the raw data to view if the data collection was successful.

clear all;close all; clc;
load('raw_data/walk_Pedri_05.mat')

%% Orientaion
X=Orientation.X;
Y=Orientation.Y;
Z=Orientation.Z;
t = second(Orientation.Timestamp);

% Plot 2D Plot
figure; clf;
subplot(2,3,1)
plot(t,X)
xlabel('t')
ylabel('X')
%xlim([0 1000])
title('2D Orientation (X)')

subplot(2,3,2)
plot(t,Y)
xlabel('t')
ylabel('Y')
%xlim([0 1000])
title('2D Orientation (Y)')

subplot(2,3,3)
plot(t,Z)
xlabel('t')
ylabel('Z')
%xlim([0 1000])
title('2D Orientation (Z)')
%% Acceleration
accelX = Acceleration.X;
accelY = Acceleration.Y;
accelY = Acceleration.Z;
accelt = second(Acceleration.Timestamp)';

%%% Plot 2d Plot
subplot(2,3,4)
plot(accelt,accelX, 'b')
title('2D Acceleration (X)')
xlabel('t')
ylabel('X')

subplot(2,3,5)
plot(accelt,accelY, 'b')
title('2D Acceleration (Y)')
xlabel('t')
ylabel('Y')

subplot(2,3,6)
plot(accelt,accelY, 'b')
title('2D Acceleration (Z)')
xlabel('t')
ylabel('Z')

