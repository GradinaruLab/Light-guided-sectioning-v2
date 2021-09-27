%only once
%run('E:\MatlabFiles\vlfeat-0.9.21\toolbox\vl_setup')
% Alon Greebaum, Cell Reports 2021

%Define parameters
close all;
path = 'E:\ImplantAnalysisData\Raw\160\';
distanceFromSurface = 516; %um
trim = 50; %um

%Go to the folder and align the images
curPath = pwd;
cd(path);
cd '.\LowExposure';
listLow = dir('*.JPG');
img1Low = imread(listLow(end).name);
img1Low = single(rgb2gray(img1Low));
cd ..
cd '.\HighExposure';
listHigh = dir('*.JPG');

%Load the last image and align all the images to it
img1 = imread(listHigh(end).name);
img1 = single(rgb2gray(img1));
figure; imshow(img1,[]);

%Pinpoint the light spot
[x, y] = getpts;
alignedMeasurmentStruct = struct('img', img1, 'imgRegHighExp', img1, 'imgRegLowExp', img1Low, 'transform', 0, 'distanceFromSurface', distanceFromSurface, 'lightCenterCord', [x; y], 'radiusLightPixels', 60);
name = [num2str(distanceFromSurface),'.mat'];
cmd = ['save ', name,' alignedMeasurmentStruct;'];
cd ..
eval(cmd);
cd '.\HighExposure';

%Load the rest and align
for ii = 1:(numel(listHigh)-1)
    img2 = imread(listHigh(ii).name);
    img2 = single(rgb2gray(img2));
    if (mod(ii,5) == 0)
        showImages = true;
    else 
        showImages = false;
    end
    %Align the images
    [img2registered, tform] = regImagesUsingSIFT(img1, img2, showImages);
    %Save to a structure
    alignedMeasurmentStruct.imgRegHighExp = img2registered;
    alignedMeasurmentStruct.transform = tform;
    alignedMeasurmentStruct.distanceFromSurface = distanceFromSurface + (numel(listHigh) - ii)*trim;
    img3 = imread(['..\LowExposure\',listLow(ii).name]);
    img3 = single(rgb2gray(img3));
    img3registered = imwarp(img3,tform,'OutputView',imref2d(size(img1)));
    alignedMeasurmentStruct.imgRegLowExp = img3registered;
    name = [num2str(alignedMeasurmentStruct.distanceFromSurface),'.mat'];
    cd ..
    cmd = ['save ', name,' alignedMeasurmentStruct;'];
    eval(cmd);
    cd '.\HighExposure';
    if (showImages)
       DrawCircle(img2registered, alignedMeasurmentStruct.lightCenterCord(1), alignedMeasurmentStruct.lightCenterCord(2), alignedMeasurmentStruct.radiusLightPixels);
       DrawCircle(img3registered, alignedMeasurmentStruct.lightCenterCord(1), alignedMeasurmentStruct.lightCenterCord(2), alignedMeasurmentStruct.radiusLightPixels);
    end
    
end

