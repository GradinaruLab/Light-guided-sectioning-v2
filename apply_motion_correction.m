function apply_motion_correction(method)
% Anat Kahan, Ryan Cho, Cell Reports 2021
% after saving all 2P file as tif using 'Zstack_show_movies_persection'
% now they are motion corrected before continue to CalmAn-Matlab
% 'demo-script' 
% run this function from the folder where the files are
clear M2 Y

listOfFiles = dir('*_v2.tif');
%Sort by date
[~, ind] = sort([listOfFiles.datenum]);
listOfFiles = listOfFiles(ind);

for fi=1:length(listOfFiles)
    Y=read_file(listOfFiles(fi).name);
    switch method
        case 'rigid'
             M2=rigid_MC(Y);
             figure;imshow(uint16(mean(M2,3)));
             imwrite(uint16(mean(M2,3)),[listOfFiles(fi).name(1:end-4) '_MC_mean.tif'])
             imwrite(uint16(max(M2,[],3)),[listOfFiles(fi).name(1:end-4) '_MC_max.tif'])
             saveastiff(M2,[listOfFiles(fi).name(1:end-4) '_v2MC.tif'])
        case 'non-rigid'
             M2=nonrigid_MC(Y);
      %       figure;imshow(uint16(mean(M2,3)));
            imwrite(uint16(mean(M2,3)),[listOfFiles(fi).name(1:end-4) '_nrMC_mean.tif'])
             imwrite(uint16(max(M2,[],3)),[listOfFiles(fi).name(1:end-4) '_nrMC_max.tif'])
            saveastiff(M2,[listOfFiles(fi).name(1:end-4) '_v2nrMC.tif'])
    end
    
end