%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Inverse model for large memory calculation
% Author: Waleed Tahir, Hao Wang
% Email: waleedt@bu.edu, wanghao6@bu.edu
% Date: 12/24/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
addpath(genpath('fcns'));   % add correponding functions

%% set parameters
spt = 25;     % z-axis sampling distance is lamdba/(4*n0) times 25, which is around 169 slices
density = 1.6;
nr = 1.59;
dn = nr-1.33;
div = 4/spt;
res = mod(4210,spt)/spt;     % 4210 is the number of layers with dz=lambda/(n0*4)
data_num = 1;
tau = 5;         % regularization term parameter
isSim = true;     % reconstruct simulated or experimental data
%% load ground truth object
if isSim==true
    % for the simulated data
    load(strcat('..\object\simulatedDownsampledData\density_', num2str(density), '.mat'));
    [nx, ny, nz] = size(data);
    f_g = data*dn;  
    f_g = gpuArray(f_g);
else
    % reconstruct object for the experimental data
    f_g = [];
end

%% load holograms
if isSim==true
    %load generated hologram data which is created by the 16840 slices, used as the measurement.
    holo_name = strcat('..\holograms\simulatedHologram\sampleHolograms\density', num2str(density), '_dn', num2str(dn),'.mat');
    load(holo_name);
    meas_g = abs(meas_g).^2;     %Note: the simulated result is complex field
else
    % load experimental data
    save_name = strcat('..\holograms\experimentalHologram\density', num2str(density),'.mat');
    load(save_name);
    meas_g = meas;
    [nx, ny] = size(meas_g);
    nz = 169;
end

%% Simulate optical setup
% create object for imaging setup
bpm = BPM_GPU(nx, ny, div);
bpm = bpm.set_nr(nr);
bpm = bpm.set_Nz(nz);
bpm = bpm.set_adjdz(res*bpm.dz);   % back-propagate for the last slice
opt = GradientOptimizer_GPU(bpm);  % create object for optimization, inverse model

%% reconstruct object

% Select reconstruction method:
r_method = 'Gradient_L1_Prox';  
switch r_method
    case 'Gradient_L1_Prox'
        tau = tau;
end

% function(measurement, ground truth of object, regularization term, density_num, d_n)
% output: recosntructed object, saving the object and images(uin8, Binary, cost)
fhat = opt.fista_tv_saving(meas_g, f_g, r_method, tau, density, bpm.dn);   
% add prior knowledge values larger than 0
mask = (fhat>0);
fhat = mask.*fhat;

%% save result
dir = sprintf('../results/Rg%1.2f/data%d/tau%02.5e/',density, dn, tau);
if ~exist(dir, 'dir')
    mkdir(dir)
end
filename = sprintf('%s/fhat.mat', dir);
save(filename,'fhat', '-v7.3');
disp('finished reconstructed one hologram');


rmpath(genpath('fcns'));   % delete added path 


