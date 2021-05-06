%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Contains functions to calculate gradient, gradient descent
%              based inverse model, showing result.
% Author: Waleed Tahir, Hao Wang
% Email: waleedt@bu.edu, wanghao6@bu.edu
% Date: 21/08/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef GradientOptimizer_GPU
    
    properties
        
        % scattering model object
        scatter_mdl;
        
        % ground truth object
        f_gt;
        
        % initialized object estimate
        f_init;
        
        % gradient descent step size
        step_size;
        
        % total number of iterations to run
        num_iter;
        
        % plot or not plot diagnostics during iterations
        is_plot_diagnostics;
        
        % save data every n iters
        save_every_nIters;
        
        % plot or not object maxprojection during reconstruction
        is_plot_object_maxproj
    end
    
    methods
        %%%%%%%%%%%%%%%% Constructor
        function this = GradientOptimizer_GPU(scatter_mdl)
            
            % scattering model object
            this.scatter_mdl = scatter_mdl;
            
            % gradient descent step size
            this.step_size = scatter_mdl.step_size;
            
            % total number of iterations to run
            this.num_iter = scatter_mdl.num_iter;
            
            % plot or not plot diagnostics during iterations
            this.is_plot_diagnostics = scatter_mdl.is_plot_diagnostics;
            
            % save data every n iters
            this.save_every_nIters = scatter_mdl.save_every_nIters;
            
            % plot or not object maxprojectin during reconstruction
            this.is_plot_object_maxproj = scatter_mdl.is_plot_object_maxproj;
            
        end

        %%%%%%%%%%%%%%%% FISTA with TV - saving results
        function [fhat, diagnostics_fig] = fista_tv_saving(this, meas_g, f_gt_g, choice, tau, density_n,d_n)
        % function(measurement, ground truth of object, regularization term, density_n, data_n)
        
        % Initialization for reconstruction   
            tvCost =  @(x) 0;
            % gradient descent diagnostics
            costs = zeros(3,this.num_iter);
            snrs = zeros(this.num_iter,1);
            step = this.step_size;
            fid_cost_prev = 0;
            % initialize object
            switch this.scatter_mdl.init
                case 'zeros'   % simple choice for initilization
                    fhat_g = zeros(this.scatter_mdl.Nx,this.scatter_mdl.Ny,this.scatter_mdl.Nz,'gpuArray');
                otherwise
                    disp('Unknown initialization');
            end              
            s_g = fhat_g;
            
        % FISTA method params
            q = 1;
            
        % Iterative reconstruct object
            for iIter = 1:this.num_iter
                tic
                % calculate gradient
                [grad_g, meas_hat_g] = this.scatter_mdl.gradient(s_g, meas_g);
                
                switch choice
                    case 'Gradient_L1_Prox'
                        % Gradient descent method + L1 proximal
                        s_g = proxL1_gpu(s_g-step*grad_g, step*tau);
                        fhatnext_g = s_g;
                        qnext = 0.5*(1+sqrt(1+4*q^2));
                        s_g = fhatnext_g + ((q-1)/qnext)*(fhatnext_g-fhat_g);        
                        q = qnext; 
                end
                
        % compute diagnostics, including different kinds of errors
                iter_time = toc;
                if isempty(f_gt_g)
                    [costs(:,iIter), ~] = this.compute_print_diagnostics(meas_g, meas_hat_g, f_gt_g, fhat_g, fhatnext_g, iIter, iter_time, tau*tvCost(fhatnext_g));
                else
                    [costs(:,iIter), snrs(iIter)] = this.compute_print_diagnostics(meas_g, meas_hat_g, f_gt_g, fhat_g, fhatnext_g, iIter, iter_time, tau*tvCost(fhatnext_g));
                end
                fid_cost_now = costs(1,1);
                
                % plot diagnostics
                % plot errors
                if this.is_plot_diagnostics
                    diagnostics_fig = this.plot_diagnostics(iIter, costs, snrs, f_gt_g);
                else
                    diagnostics_fig = [];
                end
                
                % plot projections of different plane
                if this.is_plot_object_maxproj
                    this.plot_object_maxproj(fhatnext_g, iIter);
                end
                
                % save results, name includes num of particles and dn
                if rem(iIter,this.save_every_nIters)== 0
                    fhatnext = gather(fhatnext_g);
                    save_results(fhatnext, diagnostics_fig, iIter, tau, density_n, d_n);
                end
                
                % update step size, object, and judge for iteration
                tol = 1e-1;
                if fid_cost_prev - fid_cost_now < tol
                   step = step / 1;
                end
                fhat_g = fhatnext_g;
                fid_cost_prev = fid_cost_now;
                         
            end
            
            fhat = gather(fhat_g);
        end
        
        
        %% help function
        %Print diagnostics, including iter num, errors for fidelty, tv-cost, total cost, snr, difference between new & old object, one iter time
        function [costs, snr] = compute_print_diagnostics(~, meas, meas_hat, f_gt, fhat, fhatnext, iIter, iter_time, tv_cost_g)
            
            tv_cost = gather(tv_cost_g);   % already times regularization params tau 
            fid_cost = gather(norm(meas_hat(:)- meas(:))^2);
            total_cost = fid_cost + tv_cost;
            costs(1,1) = fid_cost;
            costs(2,1) = tv_cost;
            costs(3,1) = total_cost;
            diff = gather(norm(fhatnext(:)-fhat(:))/norm(fhat(:)));
            
            if isempty(f_gt)
                snr = [];
                fprintf('[iter = %d] [fid = %e] [tv = %e] [tcost = %e] [diff. = %8.4f] [time = %f]\n', iIter, fid_cost, tv_cost, total_cost, diff, iter_time);
            else
                snr = gather(20*log10(norm(f_gt(:))/norm(f_gt(:)-fhatnext(:))));
                fprintf('[iter = %d] [fid = %e] [tv = %e] [tcost = %e] [snr = %3.4f dB] [diff. = %8.4f] [time = %f]\n', iIter, fid_cost, tv_cost, total_cost, snr, diff, iter_time);
            end

        end
        
        %%%%%%%%%%%%%%%% Plot diagnostics for errors/ snr for each step
        function fig_handle = plot_diagnostics(this, iIter, costs, snrs, f_gt_g)
            fig_handle = figure(223);
            set(223, 'Name', sprintf('t = %d', iIter));
            if isempty(f_gt_g)
                semilogy(1:iIter, costs(1,1:iIter), 'b-');
                hold on;
                o1 = semilogy(iIter, costs(1,iIter), 'bo');
                semilogy(1:iIter, costs(2,1:iIter), 'r-');
                o2 = semilogy(iIter, costs(2,iIter), 'ro');
                semilogy(1:iIter, costs(3,1:iIter), 'g-');
                o3 = semilogy(iIter, costs(3,iIter), 'go');
                hold off
                set(get(get(o1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                set(get(get(o2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                set(get(get(o3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                legend('fid cost','tv cost','total cost');
                xlim([1 this.num_iter]);
                xlabel('cost');
                ylabel('Iteration #')
                grid on;
                title(sprintf('Fid Cost: %.2e', costs(1,iIter)));
            else
                subplot(1, 3, 1:2);
                semilogy(1:iIter, costs(1,1:iIter), 'b-',...
                    iIter, costs(1,iIter), 'ro');
                xlim([1 this.num_iter]);
                grid on;
                title(sprintf('Fidelity Cost: %.2e', costs(1,iIter)));
                
                subplot(1, 3, 3);
                plot(1:iIter, snrs(1:iIter), 'b-',...
                    iIter, snrs(iIter), 'ro');
                xlim([1 this.num_iter]);
                grid on;
                title(sprintf('SNR: %.2f dB', snrs(iIter)));
            end
            
            drawnow;
        end
        
        %%%%%%%%%%%%%%%%%% Plot the reconstruct object projection for each step
        function plot_object_maxproj(this, fhatnext, iIter)
            fig_handle = figure(224);
            set(224, 'Name', sprintf('t = %d', iIter));
            subplot(1, 3, 1); % z-projection
            caxis([0,1e-5]);
            imagesc(squeeze(mean(fhatnext, 3)));
            xlabel('x'); ylabel('y'); title('Pz'); axis image;
            subplot(1, 3, 2); % x-projection
            caxis([0,1e-5]);
            imagesc(squeeze(mean(fhatnext, 2)));
            xlabel('z'); ylabel('y'); title('Px'); %axis equal;
            subplot(1, 3, 3); % y-projection
            imagesc(squeeze(mean(fhatnext, 1)));
            xlabel('z'); ylabel('x'); title('Py'); %axis equal;
            colorbar
        end
    end
end
