%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: User parameters
% TODO: Edit parameters to match your system
% Author: Waleed Tahir
% Email: waleedt@bu.edu
% Date: 09 Apr 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef SystemParams

    properties(Constant)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Imaging setup
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lambda0 = 632e-9;  % free space wavelength (m)
        mag = 20;          % magnification
        px_size = 3.45e-6; % camera pixel size (m)
        NA = 0.405;          % NA of objective lens
        is_pupil = 1;      % does the system have a pupil function (overriddes NA)
        n0 = 1.33;         % background refractive index
        z_cam = 2e-6;      % dist from camera (m)
        in_field_amp = 1;  % amplitude of the incident field
        is_rm_alias = 1;   % perform or not to perform padding to remove 
                           % aliasing due to circular-conv artifacts        
        is_crop_each = 0;  % set as 0, don't use crop for each.
                           % 1 for crop for each step, 0 for crop only once                   

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fwd model Parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        meas_type = 'intensity'; % measurement type: 'amplitude', 
                                 % 'intensity', 'complex'
                                 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FISTA
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        step_size = 2e-6;
        num_iter = 500;
        init = 'zeros';   % initialze with zeros (can be changed for initilizing with better object function)
        is_plot_diagnostics = 1;
        is_plot_object_maxproj = 1; % only works for fista_tv
        save_every_nIters = 100;
        profiling = false;
    end
    
    properties 
        Nx;
        Ny;
        Nz;
        dz;
        nr;
        pad_len;
        z_adjust;
    end
    
    methods
        %%%%%%%%%%%%%%%% Constructor
        function this = SystemParams()
            if this.is_rm_alias==0
                this.pad_len = 1;
            else
                if this.is_crop_each
                    this.pad_len = 1.5;
                else
                    this.pad_len = 2;
                end
            end
        end
    end
end
