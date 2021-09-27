% Anat Kahan, Cell Reports 2021
all_mice={'WT58N_LGS','WT35L_LGS','WT36R_LGS','WT316RR_LGS','WT_242_LGS','Drd1a_1816L_LGS'};%,'WT2170_LGS','Drd1_1N_LGS'
clear X_tissue X_GRIN_fixed X_GRIN_invivo  angle


comp=2
for mi=1:length(all_mice)
    [X_tissue{mi},X_GRIN_fixed{mi},X_GRIN_invivo{mi},new_scale{mi},angle{mi}]=extract_scale_LGS_registration(all_mice{mi},comp);
end


%% plot all scales, normalized to different imaging conditions 
figure
subplot(2,1,1)
for mi=1:length(all_mice)
    plot(X_tissue{mi},new_scale{mi},'-*','linewidth',4,'Color',mi*[0.1 0.1 0.1])
    hold on
end
xlabel('working distance (um)')
ylabel('fit scale')

subplot(2,1,2)
for mi=1:length(all_mice)
    plot(X_tissue{mi},1./new_scale{mi},'-*','linewidth',4,'Color',mi*[0.1 0.1 0.1])
    hold on
end

xlabel('working distance (um)')
ylabel('GRIN lens magnification')
xlim([0 200])

%% plot angle
figure
subplot(2,1,1)
for mi=1:length(all_mice)
    plot(X_tissue{mi},angle{mi},'-*','linewidth',4,'Color',mi*[0.1 0.1 0.1])
    hold on
end
xlabel('distance from GRIN lens surface (um)')
ylabel('angle')

%% plot distance
figure
subplot(2,1,1)
for mi=1:length(all_mice)
    plot(X_GRIN_fixed{mi},X_GRIN_invivo{mi},'-*','linewidth',4,'Color',mi*[0.1 0.1 0.1])
    hold on
end
plot([0 950],[0 950],'--k')
xlabel('imaging distance fixed (um)')
ylabel('imaging distance in vivo (um)')
ylim([450 950])
xlim([450 950])

subplot(2,1,2)
for mi=1:length(all_mice)
    plot(X_tissue{mi},X_GRIN_invivo{mi},'-*','linewidth',4,'Color',mi*[0.1 0.1 0.1])
    hold on
end
xlabel('working distance fixed  (um)')
ylabel('imaging distance in vivo (um)')
ylim([450 950])
xlim([0 200])

for mi=1:length(all_mice)
    [X_tissue2{mi},SSIM{mi}]=extract_Rcc_LGS_registration(all_mice{mi},comp);
end

%% plot SIMM-corr
figure

for ri=1:2
    subplot(2,1,ri)
    for mi=1:length(all_mice)
        this_ssim=SSIM{mi};
        plot(X_tissue2{mi},this_ssim(ri,:),'-*','linewidth',4,'Color',mi*[0.1 0.1 0.1])
        hold on
    end
    plot([-350 0],[0.5 0.5])
    ylim([0 1])
end
xlabel('distance from GRIN lens surface (um)')
ylabel('SIMM corr')

