%%%%%%%%%%%% Generate propagation function for last slice %%%%%%%%%%%%%%%%%
% input:
% Nx:   number of pixels along x-axis
% Ny:   number of pixels along y-axis
% Nz:   number of slices along z-axis
% lambda:   wavelengh in environment (lambda0/n0)
% deltaX:   pitch size of x-axis
% deltaY:   pitch size of y-axis
% deltaZ:   propagation distance of z-axis to last slice
% is_rm_alias:   control padding for object
% pad_len:  times of object size padding for calculation
% output:
% prop_last:    propagate forward in object to last slice
% prop_last2:   propagate backward in object from last slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [prop_last, prop_last2]= gen_tf_last(Nx,Ny,lambda,...
    deltaX,deltaY,deltaZ,is_rm_alias,pad_len)
%% Zero padding for anti aliasing
if is_rm_alias
    Nx = Nx * pad_len;
    Ny = Ny * pad_len;
end

%% Spatial frequency variables
k=1/lambda;
KX=[ceil(-Nx/2):1:ceil(Nx/2-1)]'.*(1/(Nx*deltaX));
KY=[ceil(-Ny/2):1:ceil(Ny/2-1)].*(1/(Ny*deltaY));
kx=repmat(KX,1,Ny); % repeat X, 1 row Ny col
ky=repmat(KY,Nx,1); % repeat Y, Nx row 1 col
kp=sqrt(kx.^2+ky.^2);
term=k.^2-kp.^2;    % frequency in propagation formula

% evanescent wave filter
eva = double(term >= 0);
term(term<0)=0;
%% Generate prop_tf for last slice
prop_last = exp(-1i*2*pi*deltaZ*(k - sqrt(term))).*eva;     % propagate forward TF in object
prop_last2 = exp(1i*2*pi*deltaZ*(k - sqrt(term))).*eva;     % propagate backward TF in object

prop_last = fftshift(prop_last);
prop_last2 = fftshift(prop_last2);
end


