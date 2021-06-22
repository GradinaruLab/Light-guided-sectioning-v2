% estimates the improvement my non rigid methods

mouse_ID={'316RR'; '316RR' ;'36R' ;'36R' ;'35L'; '35L'; '58N' ;'58N'};
mouse_depth={'810'; '720'; '755'; '845' ;'760' ;'805'; '760' ;'655'};

% can access ssim by open 'Delaunay_results' in the relevant folder 

affine=[ 0.595 0.631 0.5451 0.4833 0.6159 0.6274 0.3651 0.3491]';
Delaunay=[0.556 0.595 0.5322 0.4782 0.5988 0.5618 0.3262 0.3520]';
non_rigid=[0.5966 0.6287 0.5565 0.4879 0.6197 0.6330 0.3758 0.3574]';

T=table(mouse_ID,mouse_depth,affine,Delaunay,non_rigid);

T.Del_improvment=(T.Delaunay-T.affine)./T.affine;
T.NR_improvment=(T.non_rigid-T.affine)./T.affine;

figure
plot(ones(1,height(T)),100*T.NR_improvment,'o'); hold on
NR_improvement=mean(100*T.NR_improvment);
NR_improvement_sem=std(100*T.NR_improvment)/sqrt(length(T.NR_improvment));
plot(2*ones(1,height(T)),100*T.Del_improvment,'or')
ylabel('% improvment in ssim score')
xlim([0.5 2.5])
xticks([ 1 2 ])
xticklabels({'non rigid' ; 'Delaunay'})