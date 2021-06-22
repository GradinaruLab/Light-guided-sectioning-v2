function Register_tiff_different_days_2P
%% 03/10/20 aim to put different days of recording into one aligned fiture
clear listOfFiles
%this_path='D:\DATA_Glab\Light_sectioning\WT2170\all_tiff_to_align';
this_path='C:\Users\anatk\Documents\Light_sectioning\WT2170\all_tiff_to_align';
this_path='C:\Users\anatk\Documents\Light_sectioning\Str39_LGS\all_tiff_to_align_combined';
cd([this_path '\']) 
listOfFiles = dir('*v2_nrMC.tif');
%Sort by date
[~, ind] = sort([listOfFiles.datenum]);
listOfFiles = listOfFiles(ind);

%Load the files one after another and create mip images
clear meanY Y all_Y size_Y
k=0;
Y_Matrix=[];
Y_Matrix_odd=[];
Y_Matrix_even=[];
for ii = 1:numel(listOfFiles)
    if contains(listOfFiles(ii).name,'.tif')
        k=k+1;
       Y=read_file([this_path '\' listOfFiles(ii).name]); 
       meanY(:,:,k)=uint16(mean(Y,3));
       Y=uint16(Y);
       saveastiff(Y,['uint16YMat_' num2str(k) '.tif'])
       %all_Y{k}=Y;
      % Y_Matrix_odd = cat(3,Y_Matrix,Y(:,:,1:2:end));% takes just odd numbers
     %  Y_Matrix_even = cat(3,Y_Matrix,Y(:,:,2:2:end));% takes just even numbers
       Y_Matrix = cat(3,Y_Matrix,Y(:,:,:));% takes just even numbers
       size_Y(k)=size(Y,3);
    end
    clear Y
end
2
%% this registration is based on motion correction. based on the fact that all files has the same dimentions
M2=rigid_MC(Y_Matrix);
saveastiff(M2,['MC.tif'])

% non rigid MC
nrM2=nonrigid_MC(Y_Matrix);
saveastiff(nrM2,['nrMC.tif'])

%% This registration isbased on the idea that move the mean will work on all the matrix:
for yi=2:size(meanY,3)-1
    transform_name=['T' num2str(yi+1) 'to' num2str(yi)];  
    % check it tform available
    if exist([transform_name '.mat'])
        MT=open([transform_name '.mat']) ;
        mytform=MT.mytform;
    else
        %% define similarity points
        
        cpselect(meanY(:,:,yi+1)*5, meanY(:,:,yi)*5)% moving /fixed
        movingPoints=evalin('base','movingPoints');
        fixedPoints=evalin('base','fixedPoints');
        %%
        
        mytform = fitgeotrans(movingPoints, fixedPoints, 'NonreflectiveSimilarity');
    end
    %% apply the tranformation to the moving figure
    %registered = imwarp(fixedImage, mytform);
    Jregistered = imwarp(meanY(:,:,yi+1),mytform,'OutputView',imref2d(size(meanY(:,:,yi))));
    for fi=1:size(all_Y{yi},3)
       Jregistered_tif(:,:,fi) = imwarp(all_Y{yi}(:,:,fi),mytform,'OutputView',imref2d(size(meanY(:,:,yi))));
    end
    fixedImage=meanY(:,:,yi);
    
    save(transform_name,'mytform','fixedImage')
    transform_name2=['Tif_' num2str(yi+1) 'to' num2str(yi) '.tif'];  
    saveastiff(Jregistered_tif,transform_name2)
    % check
    tmp=uint16(mean(Jregistered_tif,3));
    figure,  subplot(1,2,1); imshowpair(meanY(:,:,yi)*5,Jregistered*5)
    subplot(1,2,2); imshowpair(tmp*5,meanY(:,:,yi)*5)
    
end


