%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: System parameters
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
        lambda0 = 0.632;   % free space wavelength (um)
        mag = 20;          % magnification
        px_size = 3.45;    % camera pixel size (um)
        NA = 0.405;        % NA of objective lens, Adjust for pupil size
        is_pupil = 1;      % does the system have a pupil function (overriddes NA)
        n0 = 1.33;         % background refractive index        
        z_cam = 2;         % dist from camera (um)
        in_field_amp = 1;  % amplitude of the incident field
        is_rm_alias = 1;   % perform or not to perform padding to remove 
                           % aliasing due to circular-conv artifacts
        pad_len = 2;     % 3 times of original object size                   
        is_crop_each = 0;  % 1 for crop for each step, 0 for crop only once                   
                                                       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fwd model Parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        meas_type = 'complex'; % measurement type: 'amplitude', 
                                 % 'intensity', 'complex'
    end
    
    properties  
        Nx;
        Ny;
        Nz;
        dz;
        nr;
    end
    
end



