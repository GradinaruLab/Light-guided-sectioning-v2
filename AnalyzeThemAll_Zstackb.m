

listOfFiles = dir('*.dat');
%Sort by date
[~, ind] = sort([listOfFiles.datenum]);
listOfFiles = listOfFiles(ind);

%%% make changes for this experiment!
averageWindowSize = 5;
start=[-0];
step=-5;
frames_per_step=1;

%Load the files one after another and create mip images
for ii = 1:numel(listOfFiles)
    temp = load2p(listOfFiles(ii).name);
    all_temp.dat=temp.dat;
    %The trivial things
   
      
     for i=1:floor(size(temp.dat,3)/frames_per_step)
           
          mkdir([listOfFiles(ii).name(1:end-4) '_' num2str(start(ii)+(i-1)*step), 'um']);
         cd([listOfFiles(ii).name(1:end-4) '_' num2str(start(ii)+(i-1)*step), 'um']);

         temp.dat=all_temp.dat(:,:,(i-1)*frames_per_step+1:i*frames_per_step);
        imwrite(max(temp.dat,[],3), ['MAX_',listOfFiles(ii).name(1:end-4), num2str(start(ii)+(i-1)*step), 'um.tif']);
        imwrite(min(temp.dat,[],3), ['MIN_',listOfFiles(ii).name(1:end-4),num2str(start(ii)+(i-1)*step), 'um.tif']);
        imwrite(uint16(round(mean(temp.dat,3))), ['MEAN_',listOfFiles(ii).name(1:end-4), num2str(start(ii)+(i-1)*step), 'um.tif']);
        imwrite(uint16(round(std(double(temp.dat),0,3))), ['STD_',listOfFiles(ii).name(1:end-4),num2str(start(ii)+(i-1)*step), 'um.tif']);
        averMean = movmean(double(temp.dat),averageWindowSize,3);
        df = double(diff(temp.dat,1,3))./averMean(:,:,2:end);
        df = df - min(df(:));
        df = df./max(df(:));
        save df df;
        
        %df = im2uint8(df);
        
        %implay(df);
    
     cd ..
     imwrite(max(temp.dat,[],3), ['MAX_',listOfFiles(ii).name(1:end-4),num2str(start(ii)+(i-1)*step), 'um.tif']);
     imwrite(uint16(round(mean(temp.dat,3))), ['MEAN_',listOfFiles(ii).name(1:end-4),num2str(start(ii)+(i-1)*step), 'um.tif']);
    end
end


for i = 1:size(temp.dat,3)
    imshow(imresize(temp.dat(:,:,i),2),[]);
    pause(.001);
end
