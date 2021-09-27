function mulitple_save_image_contours_to_behavior(A,mouse,sess,frames_num,exp,fname)
% Anat Kahan, Cell Reports 2021
%depth='755';
%section='26';
%path='D:\DATA_Glab\Light_sectioning\';
path='C:\Users\anatk\Documents\Light_sectioning\';
switch exp
    case 'behavior'
        this_dir='all_tiff_to_align';
    case 'behavior_com'
        this_dir='all_tiff_to_align_combined';
    case 'Zstack'
        this_dir='tiff_Zstack_to_align';
end
        

if ~exist(fname,'dir'); mkdir (fname);end
cd ([path mouse '\' this_dir '\sess' num2str(sess) '\' fname '\'])
for i=1:frames_num   
    imwrite(A,[path mouse '\' this_dir '\sess' num2str(sess) '\' fname '\' num2str(i) '.tif']);
end
cd('..') 
