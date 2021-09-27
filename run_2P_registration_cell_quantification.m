function run_2P_registration_cell_quantification
%LiGS Paper. Anat Kahan, 2021 
%get cell quantification : #cells in-vivo that were found in fixed
%versus # cells in-vivo total. based on manually selected cells, by mean
%projection images, by morphology (not by activity)
% uses ssim : Structural similarity (SSIM) index for measuring image
% quality, 'Exponents' — Exponents for luminance, contrast, and structural
% terms [0.8 1 1]
All_ID={'WT58N_LGS','WT35L_LGS','WT36R_LGS','WT316RR_LGS','Drd1a_1816L_LGS','WT_242_LGS'};
thresh=0.5;
%All_ID={'WT316RR_LGS'};
%All_ID={'WT36R_LGS'};
for mi=1:length(All_ID)
    mouse=All_ID{mi}
    switch mouse
        case 'WT58N_LGS'; num_reg_points=5;
        case 'WT35L_LGS'; num_reg_points=4;
        case 'WT36R_LGS';num_reg_points=5;
        case 'WT316RR_LGS';num_reg_points=6;
        case 'Drd1a_1816L_LGS';num_reg_points=4;
        case 'WT_242_LGS';num_reg_points=6;
    end
    
    for i=1:num_reg_points
        [image_score(i),score{i},is_cell{i},shuffeled_score{i}]=LGS_2P_registration_cell_quantification(mouse,i,thresh);
        sh_positive_cells(i)=length(find(shuffeled_score{i}(find(is_cell{i}))>thresh));
        positive_cells(i)=length(find(score{i}(find(is_cell{i}))>thresh));
        total_cells(i)=length(score{i}(find(is_cell{i})));
    end
     sh_positive_cell_my_mouse(mi)=sum(sh_positive_cells);% shuffeled cells- to check relayability of method
    positive_cell_my_mouse(mi)=sum(positive_cells);
    total_cell_by_mouse(mi)=sum(total_cells);
    mean_cells_per_section(mi)=mean(total_cells);
    per_positive{mi}=100*positive_cells./total_cells;
    [max_cells_per_section(mi),max_i]=max(total_cells);
    identified_cells_best_section(mi)=positive_cells(max_i);
    per_shuff_positive{mi}=100*sh_positive_cells./total_cells;
    disp([num2str(mean(per_positive{mi})) '+-' num2str(std(per_positive{mi})/sqrt(length(per_positive{mi}))) '%'])
    disp([num2str(mean(mean_cells_per_section(mi))) '+-' num2str(std(mean_cells_per_section(mi))/sqrt(length(mean_cells_per_section(mi)))) ' cells'])
   
end


for mi=1:length(All_ID)
    mean_per_positive(mi)=mean(per_positive{mi});
    med_per_positive(mi)=median(per_positive{mi});
    disp([All_ID{mi}  ' mean percent of cell identified: ' num2str(mean(per_positive{mi}))])
    disp([All_ID{mi}  ' median percent of cell identified: ' num2str(median(per_positive{mi}))])
end

disp (['total mean of mean' num2str(mean(mean_per_positive))])
disp (['total mean of median ' num2str(mean(med_per_positive))])
disp(['amount of identified cells at best section ' num2str(mean(identified_cells_best_section)) '+-' num2str(std(identified_cells_best_section)/sqrt(length(identified_cells_best_section)))]); 


figure
bh=bar(mean_per_positive); bh.FaceColor=[0.5 0.5 0.5]; hold on
for mi=1:length(All_ID)
    plot(mi*ones(1,length(per_positive{mi})), per_positive{mi},'ok'); hold on
end
%bar(per_shuff_positive)
title('percent of cells that were found to be the same cells')
xlabel('samples')
ylabel('% cells identified')
ylim([0 105])

figure
bh=bar(1, mean(per_positive)); bh.FaceColor=[0.5 0.5 0.5];
hold on 
plot(ones(1,length(per_positive)),per_positive,'ok')
hold on 
%bar(2, mean(per_shuff_positive))
title('percent of cells that were found to be the same cells')
xlabel('samples')
ylabel('% cells identified')
ylim([0 100])

disp(['mean number of identified cells: ' num2str(mean(positive_cell_my_mouse)) ' out of ' num2str(mean(total_cell_by_mouse))]); 
1
