y=load2p;
Sdata=size(y.dat,3);


%figure
for i=1:Sdata
    imwrite(y.dat(:,:,i), ['section_' num2str(i) '.tif']);
%     subplot(3,ceil(12/3),i)
%     j=i+36;
%     imagesc(y.dat(:,:,j))
%     colormap('gray')
%     title(['section ' num2str(j)])
end
