classdef BoundClass < InverseROMClass
    properties
       error_bound_Mean
       error_bound_Max
       actual_error
       error_bound_avg
       actual_error_avg
       actual_error_max
       bound_plot_mean
       bound_plot_max
       error_plot_mean
       error_plot_max
       max_snap
    end
   
    methods
        function obj = BoundClass(varargin)
        
           %obj = obj.processArgs(varargin);
           obj@InverseROMClass(varargin);
           
        end
        
        function obj = loadFullLF(obj)
            fprintf(['Loading full LF for pattern ' num2str(obj.injection) '...'])
            load([obj.top '/Results/ROM/other/LF_EIT_' num2str(obj.injection) '.mat'],'ROM')
            obj.LF{obj.injection}.V = ROM.V;
            obj.LF{obj.injection}.delta_Mean = ROM.delta_Mean;
            obj.LF{obj.injection}.delta_Max = ROM.delta_Max;
            fprintf('Done. \n')
            clear ROM
        end
       
        function obj = makeBound(obj)

            obj = obj.loadMeasurements(obj.injection);
           
            obj = obj.loadSinks();
           
            obj = obj.setUp();
            
            obj = obj.loadFullLF(); 
            
            obj.error_bound_Mean = real(obj.LF{obj.injection}.delta_Mean);
            obj.error_bound_Max = obj.LF{obj.injection}.delta_Max;
            
            for i = 1:obj.min_snap
                disp(['Finding solutions for ' num2str(i) ' snapshots...'])
                t1 = tic;
                obj.snap = i;
                mu = obj.synth_cond(obj.te);
                f = obj.RBerror(mu);
                disp(['Time after EIT Function ' num2str(toc(t1))])
                obj.actual_error = [obj.actual_error f];                 
                disp(['Relative Error between full and RB: ' num2str(f)])
            end
            
            obj = obj.cleanBound();
            
        end
        
        function f = RBerror(obj,cond_te)
            
	    mu_a = onj.makeMu(cond_te);

	    [zNh,zN] = obj.RBapprox(obj.injection,mu_a);

	    f=norm(obj.u{obj.injection}-zNh)/norm(zN);
        end
        
        
        function obj = cleanBound(obj)
            emean = obj.error_bound_Mean;
            emax = obj.error_bound_Max;
            eact = obj.actual_error;
            obj = struct();
            obj.error_bound_Mean = emean;
            obj.error_bound_Max = emax;
            obj.actual_error = eact;
        end
        
        function obj = collectBound(obj)
            % load bound class files
            for i=1:obj.num_samples
                OrderedModelClass.changePath(['Result' num2str(i)])
                top = getenv("ROMEG_TOP");
                for ii=1:obj.num_patterns
                    disp(['Loading errors for pattern ' num2str(ii)])
                    load([top '/Results/bound/bound_' num2str(ii) '.mat'],'bound')
                    if length(bound.error_bound_Mean) < obj.max_snap
                        obj.error_bound_Mean(ii,:) = NaN(1,obj.max_snap);
                        obj.error_bound_Max(ii,:) = NaN(1,obj.max_snap);
                        obj.actual_error(ii,:,i) = NaN(1,obj.max_snap);
                        clear bound
                    else
                        obj.error_bound_Mean(ii,:) = bound.error_bound_Mean(:,1:obj.max_snap);
                        obj.error_bound_Max(ii,:) = bound.error_bound_Max(:,1:obj.max_snap);
                        obj.actual_error(ii,:,i) = bound.actual_error(:,1:obj.max_snap);                        
                    end

                    clear bound
                end
            end 
        end
        
        function obj = processBound(obj)
            % organise
            %obj.error_bound_avg = mean(obj.error_bound,3,'omitnan');
            obj.actual_error_avg = mean(obj.actual_error,3,'omitnan');
            obj.actual_error_max = mean(obj.actual_error,1,'omitnan');

            obj.bound_plot_mean = mean(obj.error_bound_Mean,1,'omitnan');
            obj.bound_plot_max = mean(obj.error_bound_Max,1,'omitnan');
            obj.error_plot_mean = mean(obj.actual_error_avg,1,'omitnan');
            obj.error_plot_max = max(obj.actual_error_max,[],3);
        end
        
        function obj = saveBound(obj)
            disp('Saving processed bound data in ROM/Results/bound_data.mat')
            bound=obj;
            OrderedModelClass.changePath('ROM'); top = getenv("ROMEG_TOP");
            save([top '/Results/bound_data.mat'],'bound')
        end
        
        function obj = plotBound(obj)
            % plot
            figure()
            range = 1:obj.max_snap;
            
            s1 = loglog(range,obj.bound_plot_mean,'-o','color','k');
            hold on
            s2 = loglog(range,obj.bound_plot_max,'--.','color','r');
            s3 = loglog(range,obj.error_plot_mean,'-^','color','k');
            s4 = loglog(range,obj.error_plot_max,'-*','color','r');
            
            grid()
            %ylabel('RE')
            xlabel('Snapshots','FontSize',15)
            legend('Avg bound','Max bound','Avg error','Max error','Interpreter','latex','FontSize',15)
            
        end
    end
end
