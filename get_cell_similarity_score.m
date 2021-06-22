function [I2,xc,yc,score,is_cell_score]=get_cell_similarity_score(fixed_GRIN,invivo_registered,xi,yi,parameters,iteration)
% this function is used to find similarity between cells , for LGS_2P_registration_cell_quantification

r=parameters.radius;

for i=1:2
    if i==1
        I=fixed_GRIN;
    elseif i==2
        I=invivo_registered;
    end
    
    %[nx,ny,d] = size(I) ;
    %[X,Y] = meshgrid(1:ny,1:nx) ;
    
    %figure; imshow(I) ;
   % hold on
    th = linspace(0,2*pi) ;
    xc = xi+r*cos(th) ;
   yc = yi+r*sin(th) ;
    %plot(xc,yc,'r') ;
    
    % Keep only points lying inside rectangular
    I2{i} = imcrop(I,[xi-r yi-r r*2 r*2]);
    %imshow(I2{i})
    I2_part{i} = imcrop(I,[xi-r/3 yi-r/3 r*0.33*2 r*0.33*2]);
end

score=ssim(I2{1},I2{2},'Exponents',parameters.exponents_array);

clear I
I=invivo_registered;
k=0;
step_direction=[1:4];

while iteration.alow_iteration==1 && k<=iteration.limit
    %[nx,ny,d] = size(I) ;
   
    for iti=step_direction(1):step_direction(end)
        switch iti
            case 1
                I2{2} = imcrop(I,[xi-r+iteration.step yi-r+iteration.step r*2 r*2]);
            case 2
                I2{2} = imcrop(I,[xi-r+iteration.step yi-r-iteration.step r*2 r*2]);
            case 3
                I2{2} = imcrop(I,[xi-r-iteration.step yi-r+iteration.step r*2 r*2]);
            case 4
                I2{2} = imcrop(I,[xi-r-iteration.step yi-r-iteration.step r*2 r*2]);
        end
        M_size1=min(size(I2{1},1),size(I2{2},1));
        M_size2=min(size(I2{1},2),size(I2{2},2));

        iteration_score(iti)=ssim(I2{1}(1:M_size1,1:M_size2),I2{2}(1:M_size1,1:M_size2),'Exponents',parameters.exponents_array);
    end
    [max_value,max_ind]=max(iteration_score);
    if max_value>score
         k=k+1;
        iteration.alow_iteration=1;
        step_direction=iti;
        score=max_value;
    else
        iteration.alow_iteration=0;
        disp(['iteration ended after ' num2str(k) ' iterations'])
    end
end

% check if in-vivo is cell by comparing middle to mean value of the
% original  image
%is_cell_score=max(max(I2_part{i}))>mean(mean(I2{2}))+3*std(mean(I2{2}));
is_cell_score=max(max(I2_part{i}))>mean(mean(invivo_registered))+3*std(mean(invivo_registered));
is_cell_score=1;
end
    
