function M2 = rigid_MC(Y,GRID)

%% Load file
% tic; Y = read_file(name); toc; % read the file (optional, you can also pass the path in the function instead of Y)
Y = single(Y);                 % convert to single precision 

%% Perform nonrigid motion correction
if nargin == 1
    GRID = [25 25];
end
%% set parameters 
options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',50,'us_fac',50,'init_batch',size(Y,3));

%% perform motion correction
tic; [M2,shifts1,template1,options_rigid] = normcorre(Y,options_rigid); toc

%% compute metrics
nnY = quantile(Y(:),0.005);
mmY = quantile(Y(:),0.995);

[cY,mY,vY] = motion_metrics(Y,10);
[cM2,mM2,vM2] = motion_metrics(M2,10);
T = length(cY);

%% plot metrics
figure;
    ax1 = subplot(2,2,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
    ax2 = subplot(2,2,2); imagesc(mM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,2,3); plot(1:T,cY,1:T,cM2); legend('raw data','non-rigid'); ylim([0 1]); title('correlation coefficients','fontsize',14,'fontweight','bold')
    subplot(2,2,4); scatter(cY,cM2); hold on; plot([0.9*min(cY),1.05*max(cM2)],[0.9*min(cY),1.05*max(cM2)],'--r'); axis square;
        xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
    linkaxes([ax1,ax2],'xy')
%% Output
