function [X_Tissue,X_GRIN_fixed,X_GRIN_invivo,new_scale,angle]=extract_scale_LGS_registration(mouse,comp)
%% this function extract the scale parameter from the registrations
% LGS project
fig=0;
file_name='LGS registration table'; % 3 is with estrus state refered to proestrus
path='D:\Data_Glab\Light_sectioning\';
full_path=[path file_name '.xlsx'];
%[NUM,TXT,RAW]=xlsread(full_path,'VS AK');
[NUM,TXT,RAW]=xlsread(full_path);

switch mouse
    case 'WT36R_LGS'; num_reg_points=8;
    case 'WT58N_LGS'; num_reg_points=4;
    case 'Drd1_1N'; num_reg_points=4;
    case 'WT35L_LGS'; num_reg_points=4;
    case 'WT316RR_LGS';num_reg_points=6;
    case 'Drd1a_1816L_LGS';num_reg_points=4;
    case 'WT_242_LGS';num_reg_points=6;
end
this_mouse_ind=find(contains(TXT(:,1),mouse(1:end-4)));
folder=['D:\DATA_Glab\Light_sectioning\' mouse '\fit_files'];

switch comp
    case 1.1
        fit_name='invivo1p1'
        scale=NUM(this_mouse_ind+1,2)/NUM(this_mouse_ind+1,4)*NUM(this_mouse_ind+2,2)/NUM(this_mouse_ind+2,4);

    case 2
        fit_name='processed_940nm'
        scale=NUM(this_mouse_ind+1,4)/NUM(this_mouse_ind+1,6)*NUM(this_mouse_ind+2,4)/NUM(this_mouse_ind+2,6);

end
cd (folder)
listOfFiles2 = dir('*FIT_NRS.mat');
k=0;
for fi=1:length(listOfFiles2)
    this_file=listOfFiles2(fi).name;
    if contains(this_file,fit_name)
        k=k+1;
        A=open(this_file);
        FIT=A.FIT;
        angle(k)=FIT.this_angle;  
        new_scale(k)=FIT.this_scale/scale;        
        X_GRIN_invivo(k)=-1*str2num(this_file(strfind(this_file,'FIT')-4:strfind(this_file,'FIT')-2));
        tmp_arr=NUM(this_mouse_ind+4:this_mouse_ind+num_reg_points+3,2);
        z_ind=find(tmp_arr==X_GRIN_invivo(k));
        X_GRIN_fixed(k)=NUM(this_mouse_ind+3+z_ind,4);
        X_Tissue(k)=NUM(this_mouse_ind+3+z_ind,6);
    end 
end

if fig
figure; 
plot(new_scale,X_Tissue,'*-')
end
1
