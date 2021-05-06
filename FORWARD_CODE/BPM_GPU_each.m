%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Contains functions to initialize the BPM PSF, and forward 
% and inverse BPM models.
% Author: Waleed Tahir
% Email: waleedt@bu.edu
% Date: 01 May 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input for object
% (Nx, Ny, dz):   Nx, Ny: pixel number of x-y axis of object
%                 dz:     for z-axis (uint: lambda0/n0)
%
% Methods available to the user:
% % Forward computation
% [meas_g, u_cam_g] = forward_gpu_largeArray(f, tile_size)
% % Set number of slices 
% this = setNz(this, Nz)
% % Set input field of object
% this = set_u_in(this, u_in)
% % Set refractive index of object
% this = set_nr(this, nr) 
% 
% Note: '_g' indicates that the array exists on the GPU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef BPM_GPU_each < SystemParams
    % BPM_GPU_each properties can be seen as subclass of SystemParams
    
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
        prop_in2; % +ve prop
        prop_cam2;% +ve prop
        pupil;
        % crop image for each step
        u_crop;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%% Constructor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function this = BPM_GPU_each(Nx, Ny, div)
            this = this@SystemParams();
            this.dz = this.lambda0/(div*this.n0);   % set the z-axis sampling distance
            this.Nx = Nx;
            this.Ny = Ny;            
            this.k0 = 2*pi / this.lambda0;   % free space wavenumber
            
            % for removing aliasing artifacts
            if this.is_rm_alias  
                nx = this.pad_len*this.Nx;
                ny = this.pad_len*this.Ny;
            else
                nx = this.Nx;
                ny = this.Ny;
            end
            
            % computing TF for propagation within the volume and in air
            dx = this.px_size / this.mag;  % pitch size along x
            dy = this.px_size / this.mag;  % pitch size along y
            % propagation in water
%           lambda = this.lambda0/this.n0;
%             [prop_in_, prop_in2_, prop_cam_, prop_cam2_, pupil] = gen_tf(...
%                 this.Nx,this.Ny,lambda,...
%                 dx,dy,this.dz,this.z_cam,this.is_rm_alias,...
%                 this.is_pupil,this.NA,this.pad_len);
            % propagation in air
            [prop_in_, prop_in2_, prop_cam_, prop_cam2_, pupil] = gen_tf_air(...
                this.Nx,this.Ny,this.lambda0,this.n0,...
                dx,dy,this.dz,this.z_cam,this.is_rm_alias,...
                this.is_pupil,this.NA,this.pad_len);
            this.prop_in = gpuArray(prop_in_);   
            this.prop_in2 = gpuArray(prop_in2_);
            this.prop_cam = gpuArray(prop_cam_);
            this.prop_cam2 = gpuArray(prop_cam2_);
            this.pupil = gpuArray(pupil);
            
            % input field, plane wave
            u_in_ = this.in_field_amp * ones(nx, ny);
            this.u_in = gpuArray(u_in_);
            
            % crop map for each step or not
            if this.is_crop_each
                crop_x = 1 - exp(-((1:nx)/nx-5/6).^2*5000)*0.1;   % add mask to eliminate boundary artifacts
                crop_y = 1 - exp(-((1:ny)/ny-5/6).^2*5000)*0.1;
                u_crop = crop_x'.*crop_y;
            else
                u_crop = ones(nx, ny);
            end
            this.u_crop = gpuArray(u_crop);
        end
        
        %%%%%%%%%%%%% Methods for setting system parameters   %%%%%%%%%%%%%
        % set parameters which maybe changed for different instance
        function this = set_Nz(this, Nz)
            this.Nz = Nz;
        end
        function this = set_u_in(this, u_in)
            this.u_in = u_in;
        end
        function this = set_nr(this, nr)
            this.nr = nr;
            this.dn = nr - this.n0; 
        end
        

        %%%%%%%%%%%%% Forward compute on GPU for large arrays  %%%%%%%%%%%%
        function [meas_g, u_cam_g, u_g] = forward_largeArray(this, f, tile_size)
            
            % divide f into tiles/secctions, each of which will be loaded
            % to the gpu separately for bpm computation
            % input:
            % this:   imaging setup
            % f:      the whole object
            % tile_size:   slices loaded into GPU in each time
            % output:
            % meas_g:   measured image on camera plane
            % u_cam_g:  measured complex field 
            % u_g:      last field in each section
               
            % split data into small sections
            f_split = split_f_in_sections(f, tile_size); 
            
            % perfrom bpm computation on each of the above tiles serially
            u_g = this.u_in;   % input field
            for i = 1:length(f_split)
                f_g = gpuArray(f_split{i});   % read in the different sections
                u_g = this.forward_gpu_no3D_slice_to_slice(f_g, u_g);
            end 
                      
            % Propagate to the camera
            tmp_g = ifft2(fft2(u_g).*this.prop_cam.*this.pupil);
            u_cam_g = tmp_g(1:this.Nx,1:this.Ny);   % crop back the image
            
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
        
        %%%%%%%%%%%%% Helper function for 'forward_gpu_largeArray()'%%%%%%%
        function u_g = forward_gpu_no3D_slice_to_slice(this, f_g, u_g)  
            % propagate field in each section
            % input:
            % this:     simulation setup
            % f_g:      section of object      
            % u_g:      input field
            % output:
            % u_g:      calculated field
            [~,~,nz] = size(f_g);
            for ind_z = 1:nz
                u_g = ifft2(fft2(u_g).*this.prop_in); % diffraction step
                if this.is_rm_alias
                    % refraction step
                    u_g = u_g.*exp(1i*this.k0*this.dn*this.dz*padarray(f_g(:,:,ind_z), [this.Nx*(this.pad_len-1),...
                        this.Ny*(this.pad_len-1)], 0,'post')); 
                    % here, I need to crop the image %
                    u_g = (u_g-1).*this.u_crop+1; 
                    %u_g = (u_g-0).*this.u_crop+0;
                else
                    u_g = u_g.*exp(1i*this.k0*this.dn*this.dz*f_g(:,:,ind_z)); 
                end
            end
        end
        
    end
    
end





