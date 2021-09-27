% Anat Kahan, Alon Greebaum, Cell Reports 2021
listOfFiles = dir('*.dat');
%Sort by date
[~, ind] = sort([listOfFiles.datenum]);
listOfFiles = listOfFiles(ind);

averageWindowSize = 5;
start=[-994];
step=1;
frames_per_step=180;% check this number each time you run it
show_type='all'
%show_type='selected'; si=11;

%Load the files one after another and create mip images
for ii = 1:numel(listOfFiles)
    temp = load2p(listOfFiles(ii).name);
    if step==1
    frames_per_step=size(temp.dat,3);% check this number each time you run it
    end
    all_temp.dat=temp.dat;
    %The trivial things
    switch show_type
        case 'all'
            for i=1:floor(size(all_temp.dat,3)/frames_per_step)
            %    figure;
               % title(['section ' num2str(i)])
                temp.dat=all_temp.dat(:,:,(i-1)*frames_per_step+1:i*frames_per_step);
                for ri=1
                    for iii = 1:size(temp.dat,3)
                   %      title(['section ' num2str(i)])
                   %     imshow(imresize(temp.dat(:,:,iii),2),[]);
                   %     pause(.003);
                    end
                end
                a=temp.dat;
            %    saveastiff(a,['section_' num2str(i) '_' num2str(start+step*(i-1)) 'um_v4.tif'])
                saveastiff(a,['s' num2str(start+step*(i-1)) 'um_v2.tif'])

            end
        case 'selected'
            i=si;
            temp.dat=all_temp.dat(:,:,(i-1)*frames_per_step+1:i*frames_per_step);
            for ri=1:3
                for iii = 1:size(temp.dat,3)
                    imshow(imresize(temp.dat(:,:,iii),2),[]);
                    pause(.003);
                end
            end
            a=temp.dat;
        %    saveastiff(a,['section_' num2str(i) '_' num2str(start+step*(i-1)) 'nm_v2.tif'])
            saveastiff(a,['s' num2str(start+step*(i-1)) 'um_v2.tif'])
    end
end
1


