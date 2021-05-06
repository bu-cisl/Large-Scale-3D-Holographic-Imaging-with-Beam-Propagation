function [prop_in_, prop_in2_, prop_cam_, prop_cam2_,Pupil, prop_half]= gen_tf(Nx,Ny,lambda0,n0,...
    deltaX,deltaY,deltaZ,z_cam,is_rm_alias,ispupil,NA,pad_len)

%% Zero padding for anti-aliasing
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
term=k.^2-kp.^2;

% evanescent wave filter
eva = double(term >= 0);
term(term<0)=0;
%% Generate Pupil
if(ispupil)
    Pupil = double(kp <= (NA/lambda0));   % pupil in air
    Pupil = fftshift(Pupil);
else
    Pupil = ones(Nx,Ny);
end

%% Generate prop_tf inside object
prop_in_ = exp(-1i*2*pi*deltaZ*(k - sqrt(term))).*eva;
prop_in2_ = exp(1i*2*pi*deltaZ*(k - sqrt(term))).*eva;
prop_half = exp(-1i*2*pi*deltaZ/2*(k - sqrt(term))).*eva;

%% Generate prop_tf_cam for prop from last slice to camera
k0 = 1/lambda0;
term=k0.^2-kp.^2;
eva = double(term >= 0);
term(term<0)=0;

prop_cam_ = exp(-1i*2*pi*z_cam*(k0 - sqrt(term))).*eva;
prop_cam2_ = exp(1i*2*pi*z_cam*(k0 - sqrt(term))).*eva;

%% Propagation function
prop_in_ = fftshift(prop_in_);
prop_in2_ = fftshift(prop_in2_);
prop_cam_ = fftshift(prop_cam_);
prop_cam2_ = fftshift(prop_cam2_);
prop_half = fftshift(prop_half);

end



