%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: BPM method calculation for forward model.
% In this code, we calcualte large scale problems, the slices is 16840 for
% 500 um.
% Author: Hao Wang, Waleed Tahir
% Email: wanghao6@bu.edu, waleedt@bu.edu
% Date: 12/24/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% set parameters used in the forward model.
div = 16;   % represents the dz = lamdba/(n0*div)
Num = 9;   % large scale problem, need to divide data into different sections
R_diff = [1.59];   % refractive index contrast
density = [1.6];   % density of particles

%% calculation of forward model
for conc = density
%load raw data, divide into different sections and save
rawobj = '..\object\simulatedData\';
objfile = strcat(rawobj, 'density_', num2str(conc), '\');

for n = R_diff  
    % use a method to set system parameter dz 
    bpm = BPM_GPU_each(1024, 1024, div); 
    % set refractive index n
    bpm = bpm.set_nr(n);
    %%%%%%%%%%%%%%%% load differnt z-axis sampled object %%%%%%%%%%%%%%%%%%
    for i=1:Num   % Num of sections for each object
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load different section of object
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        savename = strcat(objfile, int2str(i), '.mat');
        load(savename);
        obj = double(data);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create imaging setup 
        [~,~,nz] = size(obj);
        bpm = bpm.set_Nz(nz);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Forward model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tile_size = 200;   % used because of large object
        [meas_g, ~, u_g] = bpm.forward_largeArray(obj, tile_size);
        bpm = bpm.set_u_in(u_g);

        i   % show which setion of object is finished 
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% View results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save and view the result
file_save = '..\holograms\simulatedHologram\';
savefile = strcat(file_save,'density_', num2str(conc),'_dn', num2str(bpm.dn), '.mat');
meas_g = gather(meas_g);  % output field
save(savefile, 'meas_g');
disp('computation finished for one data')
end
end

