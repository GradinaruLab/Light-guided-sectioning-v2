function extract_gaussian_fit_features_LiGS_slicing_by_diameter
% LiGS paper 2021, Anat Kahan and Alon Greenbaum
% this function look at the gaussian fits of samples that had implants with the same
% diameter 
% first run AnalyseStackOfImages
close all
all_diameters=[200 400 600];
for di=1:length(all_diameters)
    clear gaussian_params distanceFromSurface sigma2Dx sigma2Dy alpha2D theta2D
    diameter=all_diameters(di);
    cd('D:\DATA_Glab\Light_sectioning\LiGS pictures cryosectioning')
    switch diameter
        case 200
            list = {'sample1','sample6_200um', 'sample7_200um', 'sample8_200um', 'sample9_200um'};
        case 400
            %list = {'features166.mat', 'features216.mat', 'features1741.mat','features1762.mat','features1787.mat'};
            list = {'185_400um', '183_400um', '181_400um','168_400um'}; %,'160_400um'};
        case 600
            list = {'sample5_600um', 'sample4_600um', 'sample3_600um'};
    end
    
    for li=1:length(list)
        
        cd(list{li})
        list2 = dir('*.mat');
        for ii=1:length(list2)
            load(list2(ii).name)
            distanceFromSurface(ii)=imageStruct.distanceFromSurface;
            sigma2Dx(ii) = imageStruct.GauParam2D(3); % sigma^2= width
            sigma2Dy(ii) = imageStruct.GauParam2D(5);
            alpha2D(ii) = imageStruct.GauParam2D(1); % amplitude
            theta2D(ii) = imageStruct.GauParam2D(6)*180/pi;
            switch diameter
                case 400% camera was set to a different magnification 6/3.6
                    sigma2Dx(ii) = sigma2Dx(ii)*1.666;
                    sigma2Dy(ii) = sigma2Dy(ii)*1.666;
            end
        end
        [sorted, ind]=sort(distanceFromSurface);
        gaussian_params{li} = array2table([sorted; alpha2D(ind); 0.5*(sigma2Dx(ind)+sigma2Dy(ind))]','VariableNames',{'distance','amplitude','mean_sigma'} );
        cd ..
    end
    for li=1:length(list)
        all_MAX_D=max(gaussian_params{li}.distance);
        all_MAX_A=max(gaussian_params{li}.amplitude);
    end
    
    % now plot
    subplot(2,3,di)
    for li=1:length(list)
        plot(gaussian_params{li}.distance, gaussian_params{li}.amplitude,'-*'); hold on
        ylim([0 1500]); xlim([0 850]);
    end
    title([num2str(diameter) ' um, amplitude'])
    xlabel('distance (um)')
    ylabel('amplitude (a.u.)')
    
    subplot(2,3,di+3)
    for li=1:length(list)
        plot(gaussian_params{li}.distance, 2.355*gaussian_params{li}.mean_sigma,'-*'); hold on
        ylim([0 2.355*80]); xlim([0 850]);
    end
    title([num2str(diameter) ' um, FWHM'])
    %  clear gaussian_params
    xlabel('distance (um)')
    ylabel('mean FWHM')
end



1