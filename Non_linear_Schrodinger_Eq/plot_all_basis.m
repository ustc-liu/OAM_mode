load basis.mat;
for indx = 1:78
    subplot(8,10,indx)
    imagesc(squeeze(abs(basis(indx,:,:)).^2));
    axis equal
end

figure
for indx = 1:78
    subplot(8,10,indx)
    imagesc(angle(squeeze(basis(indx,:,:))));
    axis equal
end