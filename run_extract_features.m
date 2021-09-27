function run_extract_features(mypath)
% Written by Anat Kahan and Alon Greenbaum Cell Reports 2021 (Figure 1)
% extract features
%Load the images and extract features 
list = dir('*.mat');
intFactor = 10; %The difference in exposure times between low and high
%pixelSize = 24.16; %The pixel size of the camera when the magnification is 3.4- Alon's iphone
pixelSize = 10.9; %The pixel size of the camera when the magnification is 6 - Min's iphone 

if nargin==0
    mypath = '168_400um';
    mypath = 'sample6_200um';
end

for ii = 1:numel(list)
    
    cmd = ['load ', list(ii).name];
    eval(cmd);
    [imageStruct] = ExtractFeatures(alignedMeasurmentStruct, intFactor);
    cmd = ['save ', list(ii).name, ' imageStruct'];
    eval(cmd);
    
end


%% Plot

features = PlotFeatures(mypath, intFactor, pixelSize);
