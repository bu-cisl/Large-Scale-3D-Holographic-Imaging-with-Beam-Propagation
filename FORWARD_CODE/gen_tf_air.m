%%%%%%%%%%%%%%%% Generate different propagation functions %%%%%%%%%%%%%%%%%
% generate different propagation function for field propagation
% input:
% Nx:   number of pixels along x-axis
% Ny:   number of pixels along y-axis
% Nz:   number of slices along z-axis
% lambda:   wavelengh in environment (lambda/n0)
% deltaX:   pitch size of x-axis
% deltaY:   pitch size of y-axis
% deltaZ:   propagation distance of z-axis
% z_cam:    propagation distance from last slice to camera
% is_rm_alias:   control padding for object
% ispupil:  control pupil adding for propagate from object to camera
% NA:       numerical aperture of system
% pad_len:  times of object size padding for calculation
% output:
% prop_in:  propagate forward in object
% prop_in2_:   propagate backward in object
% prop_cam_:   propagate forward from last slice to camera
% prop_cam2_:  propagate backward from camera to last slice
% prop_fs_:
% Pupil:   pupil function for system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [prop_in_, prop_in2_, prop_cam_, prop_cam2_, Pupil]= gen_tf_air(Nx,Ny,lambda0,n0,...
    deltaX,deltaY,deltaZ,z_cam,is_rm_alias,ispupil,NA,pad_len)
%% Zero padding for anti aliasing
if is_rm_alias
    Nx = Nx * pad_len;
    Ny = Ny * pad_len;
end

%% Spatial frequency variables
lambda = lambda0/n0;
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

%% Generate Pupil
if(ispupil)
    Pupil = double(kp <= (NA/lambda0));
    p = Pupil;
    Pupil = fftshift(Pupil);   % fftshift is for calculation
else
    Pupil = ones(Nx,Ny);
end

%% Generate pupil mask for apodization, used for comparsion and adjustment.
% isapod = 1;
% if (isapod)
%     ap = 1-lambda0^2*kp.^2;
%     ap(ap<=0)=0;
%     ap = ap.^(1/4);
%     ap = fftshift(ap);
%     Pupil = ap.*Pupil;
% end

%% Generate pupil mask for water-glass-air slab, used for comparsion and adjustment.
% isslab = 1;
% if (isslab)
%     theta = 1-lambda0^2*kp.^2;
%     theta(theta<=0.1)=0.1;   
%     theta = 1-theta;
%     theta = theta.^(1/2);     % sin(theta)
%     %theta = theta/n0;
%     %theta = asin(theta)/pi*180;
%     [~, aslab] = adjust_for_single_slab(1.33, 1.33, 1.33, 500, 0.632/1.33, theta);
%     aslab = fftshift(aslab);
%     %Pupil = abs(aslab).*Pupil;
%     Pupil = (aslab).*Pupil;
% end

%% Generate prop_tf_inside
prop_in_ = exp(-1i*2*pi*deltaZ*(k - sqrt(term))).*eva;     % propagate forward TF in object
prop_in2_ = exp(1i*2*pi*deltaZ*(k - sqrt(term))).*eva;     % propagate backward TF in object

% prop_in_ = exp(-1i*2*pi*deltaZ*(0 - sqrt(term))).*eva;     % propagate forward TF in object
% prop_in2_ = exp(1i*2*pi*deltaZ*(0 - sqrt(term))).*eva;     % propagate backward TF in object

%% Generate prop_tf_cam in air
k0 = 1/lambda0;
term=k0.^2-kp.^2;
eva = double(term >= 0);
term(term<0)=0;

prop_cam_ = exp(-1i*2*pi*z_cam*(k0 - sqrt(term))).*eva;     % propagate forward TF to camera
prop_cam2_ = exp(1i*2*pi*z_cam*(k0 - sqrt(term))).*eva;     % propagate backward TF to camera

%% Generate prop_tf_firstslice (back to first slice)
% z_fs = z_cam + (Nz-1)*deltaZ;
% prop_fs_ = exp(1i*2*pi*z_fs*(k - sqrt(term))) .* fftshift(Pupil).*eva;   % propogate forward, it should be minus "-"

prop_in_ = fftshift(prop_in_);
prop_in2_ = fftshift(prop_in2_);
prop_cam_ = fftshift(prop_cam_);
prop_cam2_ = fftshift(prop_cam2_);
% prop_fs_ = fftshift(prop_fs_);
end


