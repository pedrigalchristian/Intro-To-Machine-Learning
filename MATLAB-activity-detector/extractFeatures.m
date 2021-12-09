%% Feature Extraction
%%% Features 1 to 3: Frequency via FFT (Orientation)
% We will obtain the frequency for each of three Orientation axes.
% Feature #1: Frequency of X-Axis for Orientation
% Feature #2: Frequency of Y-Axis for Orientation
% Feature #3: Frequency of Z-Axis for Orientation

%%% Features 4 to 6: Average Peak-to-Peak Amplitude via Range, Min, Max
% Feature #4: Average Peak-to-Peak Amplitude X-Axis for Orientation
% Feature #5: Average Peak-to-Peak Amplitude Y-Axis for Orientation
% Feature #6: Average Peak-to-Peak Amplitude Z-Axis for Orientation

%%% Features 7 to 9: Jerks from the Three Acceleration Axes
% Feature 7: Jerk from the X Axis
% Feature 8: Jerk from the Y Axis
% Feature 9: Jerk from the Z Axis

load('finalData_20211205_v01.mat');
JumpRopeOrientCell = JumpRopeOrientCell(1:height(JumpRopeAccelCell));


AngFreqX = [];
AngFreqY = [];
AngFreqZ = [];

AmpX = [];
AmpY = [];
AmpZ = [];

MedJerkX = [];
MedJerkY = [];
MedJerkZ = [];

LabelCol = [];

allOrientCells = whos('*OrientCell');
allAccelCells = whos('*AccelCell');

%%% Extract Features from Orientation Data
for j = 1:length(allOrientCells)
    ExerciseCell = eval(allOrientCells(j).name);
    Label = string(extractBefore(allOrientCells(j).name, "Orient"));
    for i = 1:length(ExerciseCell) 
        [Feature1, Feature2, Feature3] = PSD(ExerciseCell{i,1});
        AngFreqX = [AngFreqX; Feature1];
        AngFreqY = [AngFreqY; Feature2];
        AngFreqZ = [AngFreqZ; Feature3];
        
        [Feature4, Feature5, Feature6] = P2PAmplitude(ExerciseCell{i,1});
        
        AmpX = [AmpX; Feature4];
        AmpY = [AmpY; Feature5];
        AmpZ = [AmpZ; Feature6];
        
        LabelCol = [LabelCol; Label];
    end
end

%%% Extract Features from Acceleration Data
for j = 1:length(allAccelCells)
    ExerciseCell = eval(allAccelCells(j).name);
    for i = 1:length(ExerciseCell)
        [Feature7, Feature8, Feature9] =  Jerk(ExerciseCell{i,1});
        
        MedJerkX = [MedJerkX; Feature7];
        MedJerkY = [MedJerkY; Feature8];
        MedJerkZ = [MedJerkZ; Feature9];
    end
    
end
%% Create Table of Features
finalTable = table(AngFreqX, AngFreqY, AngFreqZ, ...
                   AmpX, AmpY, AmpZ,...
                   MedJerkX, MedJerkY, MedJerkZ, ...
                    LabelCol);
                
%%% Rescale each column of Table between 0 and 100

rescaleTable = finalTable;
for i = 1:(width(finalTable)-1)
    rescaleTable{:,i} = rescale(finalTable{:,i}, 0, 100);
end

%% Save to File
% Uncomment for either unscaled or (preferred) scaled feature data.

% writetable(finalTable, 'featureTable_20211205_v00.csv')
writetable(rescaleTable,'rescaled_featureTable_20211205_v00.csv')
 
% Debugging
%{
%{
%%
%{
clc;

dt = 1/50;
n = length(t);
freq = 1/(dt*n)*(0:n);
L = 1:floor(n/2);
figure(3);

% X
orientX = detrend(orientX);
fx = fft(orientX,n);
PSDx = fx.*conj(fx)/n;
[highestX,idX] = max(PSDx);
speedX = freq(idX)

subplot(2,3,1)
plot(freq(L),PSDx(L)) 
title('PSD of X, Max Freq: ',speedX)
xlabel('f (Hz)')
ylabel('PSD (dB/Hz)')
xlim([-1 3])

% Y
orientY = detrend(orientY);
fy = fft(orientY,n);
PSDy = fy.*conj(fy)/n;
[highestY,idY] = max(PSDy);
speedY = freq(idY)

subplot(2,3,2)
plot(freq(L),PSDy(L)) 
title('PSD of Y, Max Freq: ',speedY)
xlabel('f (Hz)')
ylabel('PSD (dB/Hz)')
xlim([-1 3])

% Z
orientZ = detrend(orientZ);
fz = fft(orientZ,n);
PSDz = fz.*conj(fz)/n;
[highestZ,idZ] = max(PSDz);
speedZ = freq(idZ)

subplot(2,3,3)
plot(freq(L),PSDz(L)) 
title('PSD of Z, Max Freq: ',speedZ)
xlabel('f (Hz)')
ylabel('PSD (dB/Hz)')
xlim([-1 3])
%}

%% Average Amplitude (Orientation)
clc;
load('walk_Pedri_01.mat')
orientX=Orientation.X;
orientY=Orientation.Y;
orientZ=Orientation.Z;
t=[1:length(orientZ)];
start = 600;
stop = 850;
orientX=orientX(start:stop);
orientY=orientY(start:stop);
orientZ=orientZ(start:stop);
t=[1:length(orientZ)];


%%
%{
% Normalize Data
normX = normalize(orientX);
normY = normalize(orientY);
normZ = normalize(orientZ);

normY = normY(1:length(t/2));
t = t/2;

% Plot to check how it looks
figure(4)
plot(t,normY)
xlabel('t')
ylabel('Y')
title('Normalized Orientation (Y)')

% Find the Average Amplitude
localmaxX = findpeaks(normX);
localmaxY = findpeaks(normY);
localmaxZ = findpeaks(normZ);

avgmaxX = mean(localmaxX);
avgmaxY = mean(localmaxY);
avgmaxZ = mean(localmaxZ);

localminXidx = islocalmin(normX);
localminYidx = islocalmin(normY);
localminZidx = islocalmin(normZ);

localMinsX = normX(localminXidx);
localMinsY = normY(localminYidx);
localMinsZ = normZ(localminZidx);

roundedLocalMinsX = round(localMins,1)

avgminX = mean(localMinsX);
avgminY = mean(localMinsY)
avgminZ = mean(localMinsZ);



% ampX = avgmaxX - avgminX
% ampY = avgmaxY - avgminY
% ampZ = avgmaxZ - avgminZ

% localMins = normY(localminY)
% roundedLocalMins = round(localMins,1)
% -mode(-roundedLocalMins)
% mode(roundedLocalMins)

%  peak2peak(normY)
%}
%}

% for i = 1:length(SquatAccelCell)
%     
% end

singleSquat = BikeAccelCell{30};
% singleSquat(:,2:4) = normalize(singleSquat(:,2:4));


squatX = singleSquat{:,2};
t = singleSquat{:,1};

lowpassFilt = lowpass(squatX, 0.5, 50);

figure(3); 
clf;
plot(t, squatX, '-k'); hold on;
plot(t, lowpassFilt, '-r')

[freq_vec, fft_vec] = OneSidedFFT(t, lowpassFilt - mean(lowpassFilt));

figure(4);
clf;
plot(freq_vec, fft_vec)


[medJerkSquat100X, medJerkSquat100Y, medJerkSquat100Z] = Jerk(singleSquat)
[freqAccelSquat100X, freqAccelSquat100Y, freqAccelSquat100Y] = PSD(singleSquat)                    
%}  

%% Functions

function [Feature1, Feature2, Feature3] = PSD(singleCell)
% Extract the frequency of Power Spectral Density (PSD) from a timetable of data
%
% Arguments:
%   SingleCell = Timetable variable from a single cell
% 
% Output:
%   Feature1 = maximum PSD for 1st column (X)
%   Feature2 = maximum PSD for 2nd column (Y)
%   Feature3 = maximum PSD for 3rd column (Z)

    % Obtain number of bins
    n = length(singleCell.Timestamp);
    dt = singleCell.Timestamp(2) - singleCell.Timestamp(1);
    freq = 1/(dt*n)*(0:n);
    L = 1:floor(n/2);
    
    singleCell.X = detrend(singleCell.X);
    singleCell.Y = detrend(singleCell.Y);
    singleCell.Y = detrend(singleCell.Y);
      
    
    % Filter out DC Component and perform FFT
    fx = fft(singleCell.X - mean(singleCell.X),n);
    fy = fft(singleCell.Y - mean(singleCell.Y),n);
    fz = fft(singleCell.Z - mean(singleCell.Z),n); 
    
    % Obtain PSD
    PSDx = fx.*conj(fx)/n;
    PSDy = fy.*conj(fy)/n;
    PSDz = fz.*conj(fz)/n;
    
    
    % Report Features
    [highestX,idX] = max(PSDx);
    Feature1 = freq(idX);

    [highestY,idY] = max(PSDy);
    Feature2 = freq(idY);

    [highestZ,idZ] = max(PSDz);
    Feature3 = freq(idZ);
    
end

function [Feature4, Feature5, Feature6] = P2PAmplitude(SingleCell)
% Extract the average peak-to-peak amplitude from a timetable of data.
%
% Arguments:
%   SingleCell = Timetable variable from a single cell
% 
% Output:
%   Feature4 = average peak-to-peak amplitude for 1st column (X)
%   Feature5 = average peak-to-peak amplitude for 2nd column (Y)
%   Feature6 = average peak-to-peak amplitude for 3rd column (Z)

     normalizedCell = SingleCell{:, 1:end};
    
    for i = 2:width(normalizedCell)
        localminXidx = islocalmin(normalizedCell(:,i), 'MaxNumExtrema', 1, 'MinSeparation', 0.25);
        localMinsX = normalizedCell(localminXidx,i);
        localmaxXidx = islocalmax(normalizedCell(:,i), 'MaxNumExtrema', 1, 'MinSeparation', 0.25);
        localMaxsX = normalizedCell(localmaxXidx,i);
        if isempty(localMinsX) || isempty(localMaxsX)
            continue
        end
        feat(i) = localMaxsX - localMinsX;
    end
    
    Feature4 = feat(2);
    Feature5 = feat(3);
    Feature6 = feat(4);
    
    % Debugging
    %{ 
    %%% Obtain Mode Minimum Peak
    % Get Unrounded and Rounded Local Mins
    LocalMinsIdx = islocalmin(normalizedCell); 
    LocalMins = normalizedCell.*LocalMinsIdx;
    roundedLocalMins = round(LocalMins,1);
    % Determine most frequent minimum peak
    minIdx = find(LocalMins == mode(roundedLocalMins),1);
    MinPeaks = LocalMins(minIdx)
    
    %%% Obtain Mode Maximum Peak
    LocalMaxsIdx = islocalmax(normalizedCell);
    LocalMaxs = normalizedCell.*LocalMaxsIdx;
    roundedLocalMaxs = round(LocalMaxs,1);
    maxIdx = find(roundedLocalMaxs == mode(roundedLocalMaxs),1);
    MaxPeaks = LocalMaxs(maxIdx)
    %}

end


function [Feature7, Feature8, Feature9] =  Jerk(SingleCell)
% Calculate the MedianX-Acceleration Gradient using Central Finite
% Difference Method.
%
% Arguments:
%   SingleCell = Timetable variable from a single cell
% 
% Output:
%   Feature7 = median time-derivative of acceleration for 1st column (X)
%   Feature8 = median time-derivative of acceleration for 2nd column (Y)
%   Feature9 = median time-derivative of acceleration for 3rd column (Z)

   
    Acceleration = SingleCell;

    % Boundary Conditions
    gradAccel(1, :) = Acceleration{1, :};
    gradAccel(height(Acceleration), :) = Acceleration{end, :};

    delta = 0.01;

    % Central Difference FDM for Time Derivative of Acceleration Data
    for k = 2:height(Acceleration) - 1
        gradAccel(k, :) =  (Acceleration{k + 1, :} - Acceleration{k - 1, :}) / (2*delta);
    end
   
    % Report Features
    Feature7 = median(gradAccel(2,:));
    Feature8 =  median(gradAccel(3,:));
    Feature9 =  median(gradAccel(4,:));
    
end