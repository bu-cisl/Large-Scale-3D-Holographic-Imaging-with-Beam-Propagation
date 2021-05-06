%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Save a 3D real numbered array to a tiff file
% Author: Waleed Tahir
% Email: waleedt@bu.edu
% Date: 12 Apr 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_results(fhat, diagnostics_fig, iIter, tau, density_n, d_n)
% find values larger than 0
mask = fhat>0;
fhat = fhat.*mask;

% normalize fhat to 0-255 (fhat_uint8)
[nx,ny,nz] = size(fhat);
fhat_uint8 = zeros(nx,ny,nz,'uint8');

for i = 1:nz
    slice = fhat(:,:,i);
    max_ = max(slice(:));
    min_ = min(slice(:));
    fhat_uint8(:,:,i) = uint8(((slice - min_)/(max_ - min_))*255);
end

% binarize fhat (fhat_bin)
thresh = 90;   % here need to consider which one is the better one.
fhat_bin = (fhat_uint8 >= thresh);

% create directory if it doesnt exist
dir = sprintf('../results/Rg%1.2f/data%d/tau%02.5e/iter_%d',density_n, d_n, tau, iIter);
if ~exist(dir, 'dir')
    mkdir(dir)
end

% save data
filename = sprintf('%s/fhat.tiff', dir);
write_mat_to_tif(fhat_uint8,filename);
filename = sprintf('%s/fhat_bin.tiff', dir);
write_mat_to_tif(fhat_bin,filename);
filename = sprintf('%s/fhat.mat', dir);
save(filename,'fhat', '-v7.3');

% save diagnostics
if ~isempty(diagnostics_fig)
    filename = sprintf('%s/cost_plot.png', dir);
    saveas(diagnostics_fig, filename);
end

end