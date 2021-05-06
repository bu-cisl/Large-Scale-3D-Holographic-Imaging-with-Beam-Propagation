function xhat = shrink3D_real(y, tau, ~)
%%% Shrinkage function for complex 3D-TV denoising.
%%% Orig Author: U. S. Kamilov
%%% Edited by: Waleed Tahir
%%% Added capability to separately threshold real and imag parts of the
%%% signal. The thresholding method can be selected via a flag.

% shortcut for computing norm along the fourth dimension
computeNormInZ = @(x) sqrt(sum(abs(x).^2, 4));


% compute the norms in dimension 4
norm_y = computeNormInZ(y);

%ampr = max(norm_y-tau, 0);
norm_y(norm_y<=0) = 1; % to avoid division by 0

xhat = repmat(max(norm_y-tau, 0)./norm_y, [1, 1, 1, 3]) .* y;

end