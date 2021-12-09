clear; close all; clc;

%%% GLOBAL variables 
sample_rate = 50; %Hz
second_of_file = 5; % seconds

%% For Walking

walk_files = dir('raw_data/walk*.mat');

WalkAccelCell = {};
WalkOrientCell = {};
for i = 1:length(walk_files)
    [outAccelCell, outOrientCell] = DivideIntoEqualSize(walk_files(i).name, second_of_file);
    
    WalkAccelCell = [WalkAccelCell; outAccelCell];
    WalkOrientCell = [WalkOrientCell; outOrientCell];
end

%%% Plot Orientation
%{
n = 700;
figure(7); plot(WalkOrientCell{n,1}.Timestamp, WalkOrientCell{n,1}.Y)
xlim('tight'); grid on;
%}

clear vars outAccelCell outOrientCell i walk_files
%% For Running

run_files = dir('raw_data/run*.mat');

RunAccelCell = {};
RunOrientCell = {};
for i = 1:length(run_files)
    [outAccelCell, outOrientCell] = DivideIntoEqualSize(run_files(i).name, second_of_file);

    RunAccelCell = [RunAccelCell; outAccelCell];
    RunOrientCell = [RunOrientCell; outOrientCell];   
end

clear vars outAccelCell outOrientCell i
%% For Biking

bike_files = dir('raw_data/bike*.mat');

BikeAccelCell = {};
BikeOrientCell = {};
for i = 1:length(bike_files)
    [outAccelCell, outOrientCell] = DivideIntoEqualSize(bike_files(i).name, second_of_file);

    BikeAccelCell = [BikeAccelCell; outAccelCell];
    BikeOrientCell = [BikeOrientCell; outOrientCell];   
end

%%% Plot Orientation
%{
figure(4); plot(BikeOrientCell{5,1}.Timestamp, BikeOrientCell{5,1}.Y)
deltaTimeOrient = BikeOrientCell{1,1}.Timestamp(end) - BikeOrientCell{1,1}.Timestamp(1);
deltaTimeAccel = BikeAccelCell{1,1}.Timestamp(end) - BikeAccelCell{1,1}.Timestamp(1);
title('Biking Modified');
ylabel('Orientation Y [deg]');
xlim('tight'); grid on; xlabel('Time [s]')
%}
clear vars outAccelCell outOrientCell i
%% For Squatting

squat_files = dir('raw_data/squat*.mat');

SquatAccelCell = {};
SquatOrientCell = {};
for i = 1:length(squat_files)
    [outAccelCell, outOrientCell] = DivideIntoEqualSize(squat_files(i).name, second_of_file);

    SquatAccelCell = [SquatAccelCell; outAccelCell];
    SquatOrientCell = [SquatOrientCell; outOrientCell];   
end

%%% Plot Orientation
%{
figure(5); plot(SquatOrientCell{1,1}.Timestamp, SquatOrientCell{1,1}.Y)
deltaTimeOrient = SquatOrientCell{1,1}.Timestamp(end) - SquatOrientCell{1,1}.Timestamp(1);
deltaTimeAccel = SquatAccelCell{1,1}.Timestamp(end) - SquatAccelCell{1,1}.Timestamp(1);
title('Squatting');
ylabel('Orientation Y [deg]');
xlim('tight'); grid on; xlabel('Time [s]')
%}

clear vars outAccelCell outOrientCell i
%% For Mountain Climbing

mtnclb_files = dir('raw_data/mountain*.mat');

MtnClbAccelCell = {};
MtnClbOrientCell = {};
for i = 1:length(mtnclb_files)
    [outAccelCell, outOrientCell] = DivideIntoEqualSize(mtnclb_files(i).name, second_of_file);

    MtnClbAccelCell = [MtnClbAccelCell; outAccelCell];
    MtnClbOrientCell = [MtnClbOrientCell; outOrientCell];   
end

%%% Plot Orientation
%{
figure(6); plot(MtnClbOrientCell{6,1}.Timestamp, MtnClbOrientCell{6,1}.Z)
deltaTimeOrient = MtnClbOrientCell{1,1}.Timestamp(end) - MtnClbOrientCell{1,1}.Timestamp(1)
deltaTimeAccel = MtnClbAccelCell{1,1}.Timestamp(end) - MtnClbAccelCell{1,1}.Timestamp(1)
title('Mountain Climbing');
ylabel('Orientation Y [deg]');
xlim('tight'); grid on; xlabel('Time [s]')
%}

clear vars outAccelCell outOrientCell i
%% For Stairs

stairs_files = dir('raw_data/stair*.mat');

StairsAccelCell = {};
StairsOrientCell = {};
for i = 1:length(stairs_files)
    [outAccelCell, outOrientCell] = DivideIntoEqualSize(stairs_files(i).name, 5);

    StairsAccelCell = [StairsAccelCell; outAccelCell];
    StairsOrientCell = [StairsOrientCell; outOrientCell];   
end

%%% Plot Orientation
%{
figure(7); plot(StairsOrientCell{5,1}.Timestamp, StairsOrientCell{5,1}.Y)
deltaTimeOrient = StairsOrientCell{1,1}.Timestamp(end) - StairsOrientCell{1,1}.Timestamp(1)
deltaTimeAccel = StairsAccelCell{1,1}.Timestamp(end) - StairsAccelCell{1,1}.Timestamp(1)
title('Stairs');
ylabel('Orientation Y [deg]');
xlim('tight'); grid on; xlabel('Time [s]')
%}

clear vars outAccelCell outOrientCell i
%% For Jumping Rope

jumprope_files = dir('raw_data/jumprope*.mat');

JumpRopeAccelCell = {};
JumpRopeOrientCell = {};
for i = 1:length(jumprope_files)
    [outAccelCell, outOrientCell] = DivideIntoEqualSize(jumprope_files(i).name, 5);

    JumpRopeAccelCell = [JumpRopeAccelCell; outAccelCell];
    JumpRopeOrientCell = [JumpRopeOrientCell; outOrientCell];   
end

%%% Plot Orientation
%{
figure(7); plot(JumpRopeOrientCell{5,1}.Timestamp, JumpRopeOrientCell{5,1}.Y)
deltaTimeOrient = JumpRopeOrientCell{1,1}.Timestamp(end) - JumpRopeOrientCell{1,1}.Timestamp(1)
deltaTimeAccel = JumpRopeAccelCell{1,1}.Timestamp(end) - JumpRopeAccelCell{1,1}.Timestamp(1)
title('Jumprope');
ylabel('Orientation Y [deg]');
xlim('tight'); grid on; xlabel('Time [s]')
%}

clear vars outAccelCell outOrientCell i
%% Save to File
% Name format: "finalData" + yyyymmdd + "v" + version number

save("finalData_20211205_v01.mat")
%% Script-Defined Functions

function [OutputAccelCell, OutputOrientCell] = DivideIntoEqualSize(MatFileName, SecondsOfFile)

    load(MatFileName);
    
    % Process Raw Data
    [~, ~, ~, ~, ~, t] = datevec(Orientation.Timestamp);
    Orientation = timetable2table(Orientation, 'ConvertRowTimes', true);
    Orientation.Timestamp = t;
      
    % Obtain Sample Period
    SamplePeriodOrient = round(Orientation.Timestamp(2) - Orientation.Timestamp(1),2)  
    % Obtain Number of Samples to achieve SecondsOfFile. This is the size
    % of one table per cell.
    NumSamplesOrient = length(0:SamplePeriodOrient:SecondsOfFile);

    
    % Remove First 5 Seconds and Last 10 Seconds of mat file
    begTrim = length(0:SamplePeriodOrient:5);
    endTrim = height(Orientation) - length(0:SamplePeriodOrient:10);
    Orientation = Orientation(begTrim:endTrim, :);
    
    % The Number of Files is how many new cells that are created
    % from the original mat file.
    NumFilesOrient = floor(height(Orientation) / NumSamplesOrient);
       
    % Connects Sinusoidal Orientation
    boolRowOrient = range(Orientation{:,1:end}) > 300; % Returns an Logical Array
    posBoolTableOrient = Orientation{:,1:end} > 0;
    negBoolTableOrient = Orientation{:,1:end} < 0;
    
    Orientation{:,1:end} = posBoolTableOrient.*(Orientation{:,1:end} - boolRowOrient.*180) ...
                          + negBoolTableOrient.*(Orientation{:,1:end} + boolRowOrient.*180);             
    
    OutputOrientCell = {};
    for i = 1:NumFilesOrient
        OutputOrientCell{end + 1, 1} = Orientation( (i-1)*NumSamplesOrient + 1 : ...
                                              (i-1)*NumSamplesOrient + ...
                                              NumSamplesOrient, ...
                                               :);
        OutputOrientCell{i, 1}.Timestamp = ...
              OutputOrientCell{i, 1}.Timestamp ...
            - OutputOrientCell{i, 1}.Timestamp(1);
       
    end
    
    % Process Raw Data
    [~, ~, ~, ~, ~, t] = datevec(Acceleration.Timestamp);
    Acceleration = timetable2table(Acceleration, 'ConvertRowTimes', true);
    Acceleration.Timestamp = t;
    
    % Obtain Sample Period
    SamplePeriodAccel = round(Acceleration.Timestamp(2) - Acceleration.Timestamp(1),2)  
    % Obtain Number of Samples to achieve SecondsOfFile. This is the size
    % of one table per cell.
    NumSamplesAccel = length(0:SamplePeriodAccel:SecondsOfFile);   
    
    % Remove First 5 Seconds and Last 10 Seconds of mat file
    Acceleration = Acceleration(begTrim:endTrim, :);
    
    % The Number of Files is how many new cells that are created
    % from the original mat file.
    NumFilesAccel = floor(height(Acceleration) / NumSamplesAccel);
    
    
    % Connects Sinusoidal Orientation
    boolRowAccel = range(Acceleration{:,1:end}) > 300; % Returns an Logical Array
    posBoolTableAccel = Acceleration{:,1:end} > 0;
    negBoolTableAccel = Acceleration{:,1:end} < 0;
    
    Acceleration{:,1:end} = posBoolTableAccel.*(Acceleration{:,1:end} - boolRowAccel.*180) ...
                          + negBoolTableAccel.*(Acceleration{:,1:end} + boolRowAccel.*180);
    
    
    % Adds Lowpass Filter for Frequencies under 0.5 Hz
    Acceleration{:,2:end} = lowpass(Acceleration{:,2:end}, 0.5, round((1/SamplePeriodAccel)) ) ;
    
    OutputAccelCell = {};
    for i = 1:NumFilesAccel
        OutputAccelCell{end + 1, 1} = Acceleration( (i-1)*NumSamplesAccel + 1 : ...
                                              (i-1)*NumSamplesAccel + ...
                                              NumSamplesAccel, ...
                                               :);
        OutputAccelCell{end, 1}.Timestamp = ...
            OutputAccelCell{end, 1}.Timestamp ...
            - OutputAccelCell{end, 1}.Timestamp(1);
        
        
    end    
    
end

  