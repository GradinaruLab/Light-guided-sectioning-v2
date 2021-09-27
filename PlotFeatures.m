function [features] = PlotFeatures(path, differenceInExposure, pixelSize)
% Written by Alon Greenbaum Cell Reports 2021 
    %Go to the folder order the files in the right order by distance from
    %surface
    curPath = pwd;
    sampleNum = path(end-2:end);
    cd(path);
    list = dir('*.mat');
    names = extractfield(list,'name');
    namesC = cellfun(@CutMat, names,'un',0);
    distanceFromSurface = cellfun(@str2num, namesC, 'un',1);
    [distanceFromSurface, ord] = sort(distanceFromSurface);
    list = list(ord);
        
    %Create the data structures to plot the data
    intensityH = zeros(1, numel(list));
    intensityL = zeros(1, numel(list));
    stdH = zeros(1, numel(list));
    stdL = zeros(1, numel(list));
    sigma = zeros(1, numel(list));
    alpha = zeros(1, numel(list));
    BFMean = zeros(1, numel(list));
    sigma2Dx = zeros(1, numel(list));
    sigma2Dy = zeros(1, numel(list));
    alpha2D = zeros(1, numel(list));
    theta2D = zeros(1, numel(list));
    
    figure;  
    ind = 0;
    
    %Create an array of axis to link them together
    numberOfAxis = 16;
    axess = gobjects(numberOfAxis);
    
    %Go one image at a time and link extract the data from the struct
    for ii = 1:numel(list)
        
        cmd = ['load ', list(ii).name];
        eval(cmd);
        %DrawCircle(imageStruct.imgRegHighExp, imageStruct.lightCenterCord(1), imageStruct.lightCenterCord(2), imageStruct.radiusLightPixels);
        intensityH(ii) = imageStruct.meanIHigh;
        intensityL(ii) = imageStruct.meanILow;
        stdH(ii) = imageStruct.stdIHigh;
        stdL(ii) = imageStruct.stdILow;
        sigma(ii) = imageStruct.sigma;
        alpha(ii) = imageStruct.alpha;
        sigma2Dx(ii) = imageStruct.GauParam2D(3);
        sigma2Dy(ii) = imageStruct.GauParam2D(5);
        alpha2D(ii) = imageStruct.GauParam2D(1);
        theta2D(ii) = imageStruct.GauParam2D(6)*180/pi;
        BFMean(ii) = mean2(imageStruct.BF);
        
        if ((ii <= numberOfAxis))
            ind = ind + 1;
            axess(ii) = subplot(sqrt(numberOfAxis), sqrt(numberOfAxis), ind); bar(imageStruct.meanArrayComMirrorX*pixelSize, imageStruct.meanArrayComMirrorY); 
            title([sampleNum,' Radial dist combined exposure, distance = ', num2str( imageStruct.distanceFromSurface),' um']);  
            xlabel('Ring location [um]'); ylabel('Mean value');
            hold on;
            ndFun = @(a, b, c, x) a*exp(-((x-b)/c).^2);
            plot(imageStruct.meanArrayComMirrorX*pixelSize, ndFun(imageStruct.alpha, imageStruct.beta,imageStruct.sigma, imageStruct.meanArrayComMirrorX), '-r', 'LineWidth',4);              
            hold off;
        end   
        
    end
    linkaxes(axess,'xy');
        
    axess2 = gobjects(numberOfAxis);    
    figure; 
    ind = 0; 
    %Draw the images of the analyzied boxes
    for ii = 1:numel(list)        
        cmd = ['load ', list(ii).name];
        eval(cmd);                
        if ((ii <= numberOfAxis))
            ind = ind + 1;
            axess2(ii) = subplot(sqrt(numberOfAxis), sqrt(numberOfAxis), ind); imagesc(imageStruct.imageBoundingBox);
            title([sampleNum,', distance = ', num2str( imageStruct.distanceFromSurface),' um']);           
        end           
    end
    linkaxes(axess2,'xy');
    
    %Fill the features structure
    features.distanceFromSurface = distanceFromSurface;
    features.intensityC = intensityL*differenceInExposure + intensityH;
    features.Sigma = sigma;
    features.Alpha = alpha;
    features.meanArrayComMirrorY = imageStruct.meanArrayComMirrorY;
    features.meanArrayComMirrorX = imageStruct.meanArrayComMirrorX*pixelSize;
    features.BFMean = BFMean; 
    features.sigma2Dx = sigma2Dx;
    features.sigma2Dy = sigma2Dy;
    features.theta2D = theta2D;
    features.alpha2D = alpha2D;
    
    
    %Plot the mean value in the circles across the stack
    axess2 = gobjects(3);
    figure; title(['Sample = ',sampleNum]);
    axess2(1) = subplot(3,1,1); stem(distanceFromSurface, intensityH); 
    title([sampleNum,' Mean value for high exposure']); xlabel('Distance from surface [um]'); ylabel('Mean value');  
    axess2(2) = subplot(3,1,2); stem(distanceFromSurface, intensityL); 
    title([sampleNum,' Mean value for low exposure']); xlabel('Distance from surface [um]'); ylabel('Mean value'); 
    axess2(3) = subplot(3,1,3); stem(distanceFromSurface, intensityL*differenceInExposure + intensityH); 
    title([sampleNum,' Mean value for combined exposure']); xlabel('Distance from surface [um]'); ylabel('Mean value'); 
    %linkaxes(axess2,'xy');
    
    %Plot the std values of the gaussian 
    ind = distanceFromSurface < 1500;
    figure; 
    subplot(2,1,1); scatter(distanceFromSurface(ind), sigma(ind)); 
    xlabel([sampleNum,' Distance from surface [um]']); ylabel('Sqrt(2) sigma');  
    subplot(2,1,2); scatter(distanceFromSurface(ind), alpha(ind)); 
    xlabel([sampleNum,' Distance from surface [um]']); ylabel('Alpha'); 
    
    cd(curPath);
end

