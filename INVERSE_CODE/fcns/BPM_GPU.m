%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Contains functions to initialize the BPM PSF, and forward 
% and inverse BPM models.
% Author: Waleed Tahir, Hao Wang
% Email: waleedt@bu.edu, wanghao6@bu.edu
% Date: 21/08/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Methods available to the user:
% % Method to generate object
% BPM_GPU(system_paras)
%
% % Forward computation
% [meas_g, u_cam_g, u_3D_g] = forward_gpu(f_g)
% [meas_g, u_cam_g] = forward_gpu_largeArray(f, tile_size)
% 
% % Inverse computation
% [grad_g, meas_hat_g] = gradient_gpu(fhat_g, meas_g)
%
% Note: '_g' indicates that the array exists on the GPU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%

classdef BPM_GPU < SystemParams
    
    properties 
        % refractive index contrast (to be used for binary [0 1] objects)
        dn;
        % input field
        u_in;
        % free space wavenumber
        k0;
        % propagation transfer function
        prop_in;
        prop_cam;
        prop_in2;      % inverse prop
        prop_cam2;     % inverse prop
        
        prop_last;     % used for the last slice
        prop_last2;    % used for the last slice
        
        prop_half;
        pupil;
        % crop image for each step
        u_crop;
    end
    
    methods
        %%%%%%%%%%%%%%%% Constructor
        function this = BPM_GPU(Nx, Ny, div)
            this = this@SystemParams;
            
            this.dz = this.lambda0/(div*this.n0);
            
            % computational grid size
            this.Nx = Nx;
            this.Ny = Ny;

            % free space wavenumber
            this.k0 = 2*pi / this.lambda0;
            
            % for removing aliasing artifacts
            if this.is_rm_alias  
                nx = this.pad_len*this.Nx;
                ny = this.pad_len*this.Ny;
            else
                nx = this.Nx;
                ny = this.Ny;
            end
            
            % computing TF for propagation within the volume and to the cam
            dx = this.px_size / this.mag;  % pitch size along x
            dy = this.px_size / this.mag;  % pitch size along y
            
            [prop_in_, prop_in2_, prop_cam_, prop_cam2_, pupil, prop_half] = gen_tf(...
                this.Nx,this.Ny,this.lambda0,this.n0,...
                dx,dy,this.dz,this.z_cam,this.is_rm_alias,...
                this.is_pupil,this.NA, this.pad_len);
            
            this.prop_in = gpuArray(prop_in_);
            this.prop_in2 = gpuArray(prop_in2_);
            this.prop_cam = gpuArray(prop_cam_);
            this.prop_cam2 = gpuArray(prop_cam2_);
            this.pupil = gpuArray(pupil);
            this.prop_half = gpuArray(prop_half);
            
            % input field
            u_in_ = this.in_field_amp * ones(nx, ny);
            this.u_in = gpuArray(u_in_);
            
            % crop map for each step or not/ maybe changed different for
            % forward/inverse model
            if this.is_crop_each
                crop_x = 1 - exp(-((1:nx)/nx-5/6).^2*5000)*0.1;   % add crop funciton to image
                crop_y = 1 - exp(-((1:ny)/ny-5/6).^2*5000)*0.1;
                u_crop = crop_x'.*crop_y;
            else
                u_crop = ones(nx, ny);
            end
            this.u_crop = gpuArray(u_crop);
        end
        
        % set parameters which maybe changed for different instance
        function this = set_Nz(this, Nz)
            this.Nz = Nz;
        end
        
        function this = set_nr(this, nr)
            this.nr = nr;
            this.dn = nr - this.n0; 
        end
        
        function this = set_u_in(this, u_in)
            this.u_in = u_in;
        end

        function this = set_adjdz(this, z_adjust)
            % this is determined by resample of object, the last slice
            % propagation, which can be adjusted accordingly
            this.z_adjust = z_adjust;
            [this.prop_last, this.prop_last2] = gen_tf_last(...
                this.Nx,this.Ny,(this.lambda0/this.n0),...
                (this.px_size / this.mag),(this.px_size / this.mag),this.z_adjust,this.is_rm_alias,...
                this.pad_len);                     
        end
        
        function [u_g_new] = forward_one_slice(this, u_g_old, slice)
            % propagte one slice 
            u_g = ifft2(fft2(u_g_old).*this.prop_in); % diffraction step
            if this.is_rm_alias
                u_g = u_g.*exp(1i*this.k0*this.dz*padarray(slice, [this.Nx*(this.pad_len-1),...
                    this.Ny*(this.pad_len-1)], 0,'post')); % refraction step
                if this.is_crop_each
                % crop for each step
                u_g = (u_g-1).*this.u_crop+1; 
                end
            else
                u_g = u_g.*exp(1i*this.k0*this.dz*slice); % refraction step
            end    
            u_g_new = u_g;
        end

        %%%%%%%%%%%%%%%% Forward computation on GPU
        function [meas_g, u_cam_g, u_3D_g] = forward(this, f_g)
            
            u_g = this.u_in;      
            
            %%%%% padded object and field after calculation
            if  this.is_rm_alias
                nx = this.Nx*this.pad_len;
                ny = this.Ny*this.pad_len;
            else
                nx = this.Nx;
                ny = this.Ny;
            end
           % record partial slices to reduce the memroy used for large
           % scale calculation. (Need to pay attention more.)
           [u_3D_g, u_g] = Partial_array([nx,ny], this.Nz-1, "complex", @(old, i)this.forward_one_slice(old, f_g(:,:,i)), u_g);

           %%% propagate for the last slice
           u_g = ifft2(fft2(u_g).*this.prop_last); % diffraction step
            if this.is_rm_alias
                u_g = u_g.*exp(1i*this.k0*this.z_adjust*padarray(f_g(:,:,this.Nz), [this.Nx*(this.pad_len-1),...
                    this.Ny*(this.pad_len-1)], 0,'post')); % refraction step
                if this.is_crop_each
                % crop for each step
                u_g = (u_g-1).*this.u_crop+1; 
                end
            else
                u_g = u_g.*exp(1i*this.k0*this.dz*f_g(:,:,this.Nz)); % refraction step
            end
  
            % Propagate to the camera
            tmp_g = ifft2(fft2(u_g).*this.prop_cam.*this.pupil);
            u_cam_g = tmp_g(1:this.Nx,1:this.Ny);
            
            switch this.meas_type
                case 'amplitude'
                     meas_g = abs(u_cam_g);
                case 'intensity'
                     meas_g = abs(u_cam_g).^2;
                case 'complex'
                     meas_g = u_cam_g;
                otherwise
                     disp('Unknown measurement method');
            end
        end
        
        %%%%%%%%%%%%%%%% Gradient computation on GPU
        function [grad_g, meas_hat_g] = gradient(this, fhat_g, meas_g)

            [meas_hat_g, u_cam_g, u_3D_g] = this.forward(fhat_g);
            
            switch this.meas_type
                case 'amplitude'
                    res_g = ((meas_hat_g - meas_g) ./ meas_hat_g) .* u_cam_g;
                case 'intensity'
                    res_g = (meas_hat_g - meas_g) .* u_cam_g;
                case 'complex'
                     res_g = meas_hat_g - meas_g;
                otherwise
                     disp('Unknown measurement method');
            end

            % parameters for padding and crop
            if  this.is_rm_alias
                nx = this.Nx*this.pad_len;
                ny = this.Ny*this.pad_len;
                F = @(x) fft2(x,nx,ny);
            else
                F = @(x) fft2(x,this.Nx,this.Ny);
            end
            
            res_g = ifft2(F(res_g).*this.prop_cam2.*this.pupil);
            %grad_g = zeros(this.Nx,this.Ny,this.Nz,'single','gpuArray');
            grad_g = zeros(this.Nx,this.Ny,this.Nz,'gpuArray');
            
            %%% modify for each step gradient
            s_g = this.u_crop.*ifft2(F(u_3D_g.get(this.Nz-1)).*this.prop_last);
            s_g = conj(s_g(1:this.Nx,1:this.Ny)).*res_g(1:this.Nx, 1:this.Ny)...
                  .*conj(1i*this.k0*this.z_adjust*exp(1i*this.k0*this.z_adjust*fhat_g(:,:,this.Nz)));
            grad_g(:,:,this.Nz) = 2*real(s_g);
                            
            % calculate second term
            res_g = res_g.*conj(this.u_crop.*exp(1i*this.k0*this.z_adjust*padarray(fhat_g(:,:,this.Nz),...
                    [this.Nx*(this.pad_len-1),this.Ny*(this.pad_len-1)], 0,'post')));
            res_g = ifft2(F(res_g) .* this.prop_last2);
                        
            %%%%%% calculate gradient for other slices
            for ind_z = this.Nz-1:-1:2
                % calcualte first term
                s_g = this.u_crop.*ifft2(F(u_3D_g.get(ind_z-1)).*this.prop_in);
                s_g = conj(s_g(1:this.Nx,1:this.Ny)).*res_g(1:this.Nx, 1:this.Ny)...
                    .*conj(1i*this.k0*this.dz*exp(1i*this.k0*this.dz*fhat_g(:,:,ind_z)));
                grad_g(:,:,ind_z) = 2*real(s_g);
                
                % calculate second term
                res_g = res_g.*conj(this.u_crop.*exp(1i*this.k0*this.dz*padarray(fhat_g(:,:,ind_z),...
                    [this.Nx*(this.pad_len-1),this.Ny*(this.pad_len-1)], 0,'post')));
                res_g = ifft2(F(res_g) .* this.prop_in2);
            end            
            
            %%%%% calculate gradient for first slice
            s_g = ifft2(fft2(this.u_in).*this.prop_in).*this.u_crop;
            s_g = conj(s_g(1:this.Nx,1:this.Ny)).*res_g(1:this.Nx, 1:this.Ny)...
                .*conj(1i*this.k0*this.dz*exp(1i*this.k0*this.dz*fhat_g(:,:,1)));
            grad_g(:,:,1) = 2*real(s_g);
        end
            
            
%             %%% gradient for the last slice
%             tmp_g = ifft2(F(u_3D_g.get(this.Nz-1)).*this.prop_last);
%             s_g = tmp_g(1:this.Nx,1:this.Ny);
%             s_g = conj(s_g).*res_g(1:this.Nx, 1:this.Ny).*conj(1i*this.k0*this.z_adjust*exp(1i*this.k0*this.z_adjust*fhat_g(:,:,this.Nz)));
%             grad_g(:,:,this.Nz) = 2*real(s_g);
% 
%             res_g = res_g.*conj(exp(1i*this.k0*this.z_adjust*padarray(fhat_g(:,:,this.Nz),...
%                 [this.Nx*(this.pad_len-1), this.Ny*(this.pad_len-1)], 0,'post')));
%             res_g = ifft2(F(res_g) .* this.prop_last2);  
%             
%             %%% gradient for other slices
%             for ind_z = (this.Nz-1):-1:2
%                 %tmp_g = ifft2(F(u_3D_g(:,:,ind_z-1)).*this.prop_in);
%                 tmp_g = ifft2(F(u_3D_g.get(ind_z-1)).*this.prop_in);
%                 s_g = tmp_g(1:this.Nx,1:this.Ny);
%                 s_g = conj(s_g).*res_g(1:this.Nx, 1:this.Ny).*conj(1i*this.k0*this.dz*exp(1i*this.k0*this.dz*fhat_g(:,:,ind_z)));
%                 grad_g(:,:,ind_z) = 2*real(s_g);
%                 
%                 res_g = res_g.*conj(exp(1i*this.k0*this.dz*padarray(fhat_g(:,:,ind_z),...
%                     [this.Nx*(this.pad_len-1), this.Ny*(this.pad_len-1)], 0,'post')));
%                 res_g = ifft2(F(res_g) .* this.prop_in2);
%             end
%             
%             uin_g = this.u_in;
%             s_g = ifft2(fft2(uin_g).*this.prop_in);
%             s_g = conj(s_g(1:this.Nx,1:this.Ny)).*res_g(1:this.Nx, 1:this.Ny)...
%                 .*conj(1i*this.k0*this.dz*exp(1i*this.k0*this.dz*fhat_g(:,:,1)));
%             grad_g(:,:,1) = 2*real(s_g);
%             
            
    end
        
end
    