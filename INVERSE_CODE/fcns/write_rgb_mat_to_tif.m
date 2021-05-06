% save 3D mat as tiff file wihtout compression
function write_rgb_mat_to_tiff(mat,filename)

% save angiogram segmentation as multi-image tiff file
imwrite(mat(:,:,:,1),filename,'compression','none');
[~,~,~,nz] = size(mat);
for i = 2:nz
   imwrite(mat(:,:,:,i),filename,'compression','none','writemode','append'); 
end

end