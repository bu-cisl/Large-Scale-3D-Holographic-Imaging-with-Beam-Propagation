classdef Gradient3D_gpu < handle
    %%% This class contains methods that allow the computation of a 3D
    %%% gradient by assuming circular boundary conditions for 3D images.
    %%% U. S. Kamilov
    
    properties
        sigSize; % size of the signal
        freqResponseDTD; % Frequency response of multTranspose(mult(x))
    end
    
    methods
        function this = Gradient3D_gpu(sigSize)
            %%% Constructor
            
            %%% Store the signal size
            this.sigSize = sigSize;
            
            %%% We next generate FFT-response of the gradient operator.
            %%% This is useful for TV
            
            %%% Signal dimensions
            Nx = sigSize(2);
            Ny = sigSize(1);
            Nz = sigSize(3);
            
            %%% Frequencies
            wx = 2*pi*(0:Nx-1)/Nx;
            wy = 2*pi*(0:Ny-1)/Ny;
            wz = 2*pi*(0:Nz-1)/Nz;
            
            [Wx, Wy, Wz] = meshgrid(wx, wy, wz);
            
            %%% 3D frequency response of the gradient
            tmp = 6 ...
                -exp(-1j*Wx)-exp(1j*Wx) ...
                -exp(-1j*Wy)-exp(1j*Wy) ...
                -exp(-1j*Wz)-exp(1j*Wz);
            this.freqResponseDTD = gpuArray(tmp);
            
        end
    end
    
    methods(Static)
        function z = mult(x)
            %%% Compute gradient: d = D*x
            
            z = zeros([size(x), 3], 'gpuArray'); % initialize
            
            filter1 = cat(1,0,-1,1);
            filter2 = cat(2,0,-1,1);
            filter3 = cat(3,0,-1,1);
            
            z(:,:,:,1) = imfilter(x,filter1,'circular');
            z(:,:,:,2) = imfilter(x,filter2,'circular');
            z(:,:,:,3) = imfilter(x,filter3,'circular');
            
        end
        
        function x = multTranspose(z)
            %%% Conpute adjoint of the gradient: x = D^T*d
            
            filter1 = cat(1,1,-1,0);
            filter2 = cat(2,1,-1,0);
            filter3 = cat(3,1,-1,0);
            
            fx = z(:,:,:,1);
            fy = z(:,:,:,2);
            fz = z(:,:,:,3);
            
            fx = imfilter(fx,filter1,'circular');
            fy = imfilter(fy,filter2,'circular');
            fz = imfilter(fz,filter3,'circular');
            
            x = fx+fy+fz;
        end
    end
end