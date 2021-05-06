% read tiff file to mat
function y = read_tif_to_mat(filename)

nz = size(imfinfo(filename),1);

for i = 1:nz
   y(:,:,i) = imread(filename,i); 
end

end
