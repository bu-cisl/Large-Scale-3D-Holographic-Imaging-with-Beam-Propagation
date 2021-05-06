%xhatnext             = proxTV(xhatnext, step*tau , gradObj);
function [xhat_g, outs] = proxTV_gpu(y_g       , lam      , gradObj)

%%% TV denoising of a complex image 
%%%
%%% U. S. Kamilov

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful shortcuts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

computeSnr = @(x, xhat) 20*log10(norm(MyV2C(x))/norm(MyV2C(x)-MyV2C(xhat)));
tvCost = @(x) sum(sum(sum(sqrt(sum(abs(gradObj.mult(x)).^2, 4)))));
computeCost = @(x) 0.5*(norm(MyV2C(x)-MyV2C(y))^2)+lam*tvCost(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numIter = 10; % Number of iterations
plotRecon = false; % Plot the progress
tol = 1e-6; % Tolerance for relative change
verbose = false; % Display messages

rho = 1; % quadratic parameter
xhat0_g = y_g; % initial image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial values
xhat_g = xhat0_g;
d_g = gradObj.mult(xhat_g);
xhatmult_g = d_g;
s_g = zeros(size(d_g), 'gpuArray');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterative reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iIter = 1:numIter
    
    % message to command line
    if(verbose)
        fprintf('[proxTV: %d/%d]', iIter, numIter);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update estimates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xhatprev_g = xhat_g; % store the current estimate before update
    
    % update dw
    d_g = shrink3D_real(xhatmult_g-s_g/rho, lam/rho, 0);

    % update xhat
    dat = y_g + rho*gradObj.multTranspose(d_g + s_g/rho);
    xhat_g = real(ifftn(fftn(dat) ./ (1 + rho*gradObj.freqResponseDTD)));
    xhatmult_g = gradObj.mult(xhat_g);
    % update s
    s_g = s_g + rho*(d_g - xhatmult_g);
    
    % check tolerance
    if(norm(xhat_g(:)-xhatprev_g(:))/norm(xhatprev_g(:)) < tol)
        break;
    end
end