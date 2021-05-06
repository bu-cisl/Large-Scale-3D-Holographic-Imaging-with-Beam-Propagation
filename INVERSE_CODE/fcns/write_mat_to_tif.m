% save 3D mat as tiff file wihtout compression
function write_mat_to_tif(mat,filename)

% save angiogram segmentation as multi-image tiff file
imwrite(mat(:,:,1),filename,'compression','none');
[~,~,nz] = size(mat);
for i = 2:nz
   imwrite(mat(:,:,i),filename,'compression','none','writemode','append'); 
end

end