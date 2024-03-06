classdef SnapShotClass < PlottingClass
    properties
        range
        TRAD_history
        ROM_history
        errors_TRAD
        errors_ROM
        errors_TRAD_split
        errors_ROM_split
        errors_TRAD_max
        errors_ROM_max
        errors_TRAD_split_max
        errors_ROM_split_max
        est_vals_avg = []
        est_vals_ROM_avg
        split_conds
        omit
        split_elec
        split_layers
        removed
        idx
        layer_names
        ROM
        lgd
    end
    
    methods
        
        function obj = SnapShotClass(varargin)
        %
        %   SnapShotClass(name1,value1,name2,value2,...)
        %
        % Arguments:
        %   range           - (int) visualize how many snapshots
        %   num_samples     - array of sample numbers to read from
        %   sample_num      - sample or samples to visualize
        %   split_conds     - don't average across all samples
        %   split_layers    - split the layers being displayed
        %   omit            - array for each sample, omit convergences that
        %                     take longer than this
        %   split_elec      - don't average across electrodes
        %   tissue          - which tissue to plot
        %   layer_names     - names of tissue layers in cell array
        %   ROM             - Only load ROM values
        %
        %
            obj = obj.processArgs(varargin);
            
            if ~isempty(obj.ROM)
                obj.split_layers = true;
            end
        
        end
        
        
        function obj = readResults(obj,c,c2)
            count=1; folder = ['Result' num2str(count)];
            OrderedModelClass.changePath(folder); obj.top = getenv("ROMEG_TOP");
            obj.results_ROM = {}; obj.results_TRAD = {};
            folder2 = 'inverse_';            
            for l = 1:size(c,2)
                folder2 = [folder2 num2str(c(1,l))];
            end
            folder3 = 'inverse_';            
            for l = 1:size(c2,2)
                folder3 = [folder3 num2str(c(1,l))];
            end
            if isempty(obj.num_samples) && isempty(obj.sample_num)
                while isfolder([obj.top '/../' folder])
                    OrderedModelClass.changePath(folder); obj.top = getenv("ROMEG_TOP");
                    load([obj.top '/Results/inverse/ROM/' folder2 '/estimate_snaps.mat'],'estimates')
                    obj.results_ROM(:,:,count) = estimates;
                    load([obj.top '/Results/inverse/TRAD/' folder3 '/estimate.mat'],'estimates','histories')
                    obj.results_TRAD(:,:,count) = histories;
                    load([obj.top '/Results/measurements/prep.mat'],'Data')
                    %layers=find(logical(Data.mu_max-Data.mu_min));
                    layers = c;
                    obj.conds(:,:,count) = Data.synth_cond(layers);
                    count=count+1;
                    folder = ['Result' num2str(count)];
                end
            else
                if ~isempty(obj.sample_num)
                    samples = obj.sample_num;
                elseif ~isempty(obj.num_samples)
                    samples = 1:obj.num_samples;
                end
                for i = 1:length(samples)
                    folder = ['Result' num2str(samples(i))];
                    OrderedModelClass.changePath(folder); obj.top = getenv("ROMEG_TOP");
                    load([obj.top '/Results/inverse/ROM/' folder2 '/estimate_snaps.mat'],'estimates')
                    obj.results_ROM(:,:,i) = estimates;
                    load([obj.top '/Results/inverse/TRAD/' folder3 '/estimate.mat'],'estimates','histories')
                    obj.results_TRAD(:,:,i) = histories;
                    load([obj.top '/Results/measurements/prep.mat'],'Data')
                    %layers=find(logical(Data.mu_max-Data.mu_min));
                    layers = c;
                    obj.conds(:,:,i) = Data.synth_cond(layers);
                end
            end
        end
        
        function obj = readROMResults(obj,c)
            count=1; folder = ['Result' num2str(count)];
            OrderedModelClass.changePath(folder); obj.top = getenv("ROMEG_TOP");
            obj.results_ROM = {};
            folder2 = 'inverse_';            
            for l = 1:size(c,2)
                folder2 = [folder2 num2str(c(1,l))];
            end
            if ~isempty(obj.sample_num)
                samples = obj.sample_num;
            elseif ~isempty(obj.num_samples)
                samples = 1:obj.num_samples;
            end
            for i = 1:length(samples)
                folder = ['Result' num2str(samples(i))];
                OrderedModelClass.changePath(folder); obj.top = getenv("ROMEG_TOP");
                load([obj.top '/Results/inverse/ROM/' folder2 '/estimate_snaps.mat'],'estimates')
                obj.results_ROM(:,:,i) = estimates;
                load([obj.top '/Results/measurements/prep.mat'],'Data')
                layers = c;
                obj.conds(:,:,i) = Data.synth_cond(layers);
            end
        end
        
        function obj = processROMResults(obj)
            obj.num_samples = size(obj.results_ROM,3);
            obj.num_patterns = size(obj.results_ROM,2);
            obj.removed = zeros(1,obj.num_samples);
            for ii = 1:obj.num_samples
                for jj = 1:obj.num_patterns
                    est_vals_ROM = obj.results_ROM(:,jj,ii);
                    est_vals_ROM = interp1(1:length(est_vals_ROM{1}),est_vals_ROM{1},1:obj.range);
                    obj.est_vals_ROM_avg(:,:,jj) = est_vals_ROM;
                    final_val = est_vals_ROM(end,:);
                    
                    logi = ~(((obj.conds(:,:,ii)-obj.conds(:,:,ii)*0.1) < final_val)==0 + ((obj.conds(:,:,ii)+obj.conds(:,:,ii)*0.1) > final_val)==0);
                    if ~logi(1)
                        obj.est_vals_ROM_avg(:,:,jj) = NaN(obj.range,size(final_val,2));
                        obj.removed(ii) = obj.removed(ii) + 1;
                        continue
                    end
                end
                obj.ROM_history(:,:,ii) = mean(obj.est_vals_ROM_avg,3,"omitnan");
                obj.ROM_RE(:,:,ii) = abs(obj.ROM_history(:,:,ii)-obj.conds(:,:,ii))./obj.conds(:,:,ii);
            end
            if isempty(obj.split_conds)
                if isempty(obj.tissue)
                    obj.errors_ROM_split = mean(obj.ROM_RE,3,"omitnan");
                    obj.errors_ROM_split_max = max(obj.ROM_RE,[],3);
                    obj.errors_ROM = mean(obj.errors_ROM_split,2,"omitnan");
                    obj.errors_ROM_max = max(obj.errors_ROM_split_max,[],2);
                else
                    obj.errors_ROM_split = mean(obj.ROM_RE,3,"omitnan");
                    obj.errors_ROM_split_max = max(obj.ROM_RE,[],3);
                    obj.errors_ROM = mean(obj.errors_ROM_split(:,obj.tissue),2,"omitnan");
                    obj.errors_ROM_max = max(obj.errors_ROM_split_max(:,obj.tissue),[],2);
                end
            elseif ~isempty(obj.split_conds)
                if isempty(obj.tissue)
                    obj.errors_ROM_split = mean(obj.ROM_RE,2,"omitnan");
                else
                    obj.errors_ROM_split = mean(obj.ROM_RE(:,obj.tissue,:),2,"omitnan");
                end
            end
        end
        
        function obj = cleanProcessed(obj)
            obj.TRAD_history = [];
            obj.ROM_history = [];
            obj.errors_TRAD = [];
            obj.errors_ROM = [];
            obj.errors_TRAD_split = [];
            obj.errors_ROM_split = [];
            obj.errors_TRAD_max = [];
            obj.errors_ROM_max = [];
            obj.errors_TRAD_split_max = [];
            obj.errors_ROM_split_max = [];
            obj.est_vals_avg = [];
            obj.est_vals_ROM_avg = [];
            obj.removed = [];
            obj.TRAD_RE = [];
            obj.ROM_RE = [];
        end
            
        function obj = processResults(obj,c2)
            
            obj = obj.cleanProcessed();
            
            obj.num_samples = size(obj.results_TRAD,3);
            obj.num_patterns = size(obj.results_TRAD,2);
            obj.removed = zeros(1,obj.num_samples);
            for ii = 1:obj.num_samples
                obj.est_vals_avg = [];
                for jj = 1:obj.num_patterns
                    cell = obj.results_TRAD(:,jj,ii);
                    try
                        cond_len = size(cell{1}(1).x,2);
                        fc_vals = extractfield(cell{1},'funccount');
                        %fc_vals = fc_vals';
                        est_vals = extractfield(cell{1},'x');
                        est_vals = reshape(est_vals,[cond_len,length(fc_vals)])';
                        final_val = est_vals(end,:);

                        final_func = fc_vals(end)*4;
                        if (final_func < 8)
                            obj.est_vals_avg(:,:,jj) = NaN((obj.range*4)/4,size(final_val,2));
                            obj.removed(ii) = obj.removed(ii) + 1;
                        else
                            est_vals = interp1(fc_vals(2:end)*4,est_vals(2:end,:),4:4:obj.range*4);
                            est_vals = fillmissing(est_vals,'constant',final_val);
                            obj.est_vals_avg(:,:,jj) = est_vals;
                        end
                    catch
                        disp(['Missing traditional estimation for electrode pair ' num2str(jj) ' in sample ' num2str(ii) '.'])
                        obj.est_vals_avg(:,:,jj) = NaN((obj.range*4)/4,size(final_val,2));
                        obj.removed(ii) = obj.removed(ii) + 1;
                    end
                    
                    est_vals_ROM = obj.results_ROM(:,jj,ii);
                    est_vals_ROM = interp1(1:length(est_vals_ROM{1}),est_vals_ROM{1},1:obj.range);
                    obj.est_vals_ROM_avg(:,:,jj) = est_vals_ROM;
                    
                    logi = ~(((obj.conds(:,1:cond_len,ii)-obj.conds(:,1:cond_len,ii)*0.1) < final_val)==0 + ((obj.conds(:,1:cond_len,ii)+obj.conds(:,1:cond_len,ii)*0.1) > final_val)==0);
                    if ~logi(1)
                        obj.est_vals_avg(:,:,jj) = NaN((obj.range*4)/4,size(final_val,2));
                        obj.removed(ii) = obj.removed(ii) + 1;
                        continue
                    end
                    if isempty(obj.omit)
                        omit = final_func;
                    elseif length(obj.omit) == 1
                        omit = obj.omit;
                    else
                        omit = obj.omit(ii);
                    end
                    if final_func > omit
                        obj.est_vals_avg(:,:,jj) = NaN(obj.range*4,size(final_val,2));
                        continue
                    end
                end
                obj.TRAD_history(:,:,ii) = mean(obj.est_vals_avg,3,"omitnan");
                obj.TRAD_RE(:,:,ii) = abs(obj.TRAD_history(:,1:cond_len,ii)-obj.conds(:,1:cond_len,ii))./obj.conds(:,1:cond_len,ii);

                obj.ROM_history(:,:,ii) = mean(obj.est_vals_ROM_avg,3,"omitnan");
                obj.ROM_RE(:,:,ii) = abs(obj.ROM_history(:,:,ii)-obj.conds(:,:,ii))./obj.conds(:,:,ii);
            end
            
            if isempty(obj.split_conds)
                if isempty(obj.tissue)
                    obj.errors_TRAD_split = mean(obj.TRAD_RE,3,"omitnan");
                    obj.errors_ROM_split = mean(obj.ROM_RE,3,"omitnan");
                    obj.errors_TRAD_split_max = max(obj.TRAD_RE,[],3);
                    obj.errors_ROM_split_max = max(obj.ROM_RE,[],3);

                    obj.errors_TRAD = mean(obj.errors_TRAD_split,2,"omitnan");
                    obj.errors_ROM = mean(obj.errors_ROM_split,2,"omitnan");
                    obj.errors_TRAD_max = max(obj.errors_TRAD_split_max,[],2);
                    obj.errors_ROM_max = max(obj.errors_ROM_split_max,[],2);
                else
                    obj.errors_TRAD_split = mean(obj.TRAD_RE,3,"omitnan");
                    obj.errors_ROM_split = mean(obj.ROM_RE,3,"omitnan");
                    obj.errors_TRAD_split_max = max(obj.TRAD_RE,[],3);
                    obj.errors_ROM_split_max = max(obj.ROM_RE,[],3);
                    
                    obj.idx = ismember(obj.tissue,c2);
                    if ~isempty(find(obj.idx,1))
                        obj.errors_TRAD = mean(obj.errors_TRAD_split(:,obj.tissue(obj.idx)),2,"omitnan");
                        obj.errors_TRAD_max = max(obj.errors_TRAD_split_max(:,obj.tissue(obj.idx)),[],2);
                    end
                    
                    obj.errors_ROM = mean(obj.errors_ROM_split(:,obj.tissue),2,"omitnan");
                    obj.errors_ROM_max = max(obj.errors_ROM_split_max(:,obj.tissue),[],2);
                end
            elseif ~isempty(obj.split_conds)
                if isempty(obj.tissue)
                    obj.errors_TRAD_split = mean(obj.TRAD_RE,2,"omitnan");
                    obj.errors_ROM_split = mean(obj.ROM_RE,2,"omitnan");
                else
                    obj.idx = ismember(obj.tissue,c2);
                    if ~isempty(find(obj.idx,1))
                        obj.errors_TRAD_split = mean(obj.TRAD_RE(:,obj.tissue(obj.idx),:),2,"omitnan");
                        obj.errors_ROM_split = mean(obj.ROM_RE(:,obj.tissue(obj.idx),:),2,"omitnan");
                    end
                end
            end
        end
        
        function plotSnapshots(obj,varargin)
            
            disp(['Removed ' num2str(obj.removed) ' electrode patterns respective to each sample.' ])
            
            figure()
            if isempty(obj.split_conds)
                
                if isempty(obj.split_layers)
                    loglog(1:obj.range,obj.errors_ROM,'-^','color','k'); l1 = 'ROM - AVG RE';
                    hold on
                    loglog(1:obj.range,obj.errors_ROM_max,':x','color','r'); l2 = 'ROM - MAX RE';
                    if ~isempty(find(obj.idx,1))
                        loglog(4:4:obj.range*4,obj.errors_TRAD(:,:),'-o','color','k'); l3 = 'Traditional - AVG RE';
                        loglog(4:4:obj.range*4,obj.errors_TRAD_max(:,:),'--.','color','r'); l4 = 'Traditional - MAX RE';
                    end
                    ylabel('RE','FontSize',17)
                    xlabel('\fontname{Arial} Number of \it{n} \times \it{n}\rm linear systems to solve','interpreter','tex','FontSize',15)
                    if ~isempty(find(obj.idx,1))
                        legend(l1,l2,l3,l4,'FontSize',13);
                    else
                        legend(l1,l2);
                    end
                    grid();
                    set(gca,'FontSize',20);
                    
                elseif ~isempty(obj.split_layers)
                    icons = [{'-diamond'},{'-square'},{'-*'},{'-^'},{'-o'},{'-x'}];
                    colors = [{"#0072BD"},{'r'},{'k'},{"#D95319"},{"#EDB120"},{"#7E2F8E"}];
                    for i = obj.tissue
                        semilogy(1:obj.range,obj.errors_ROM_split(:,i),icons{i},'color',colors{i},'DisplayName',obj.layer_names{i});
                        hold on
%                         if obj.idx(i)
%                             semilogy(1:obj.range,obj.errors_TRAD_split(:,i),'-o','color','k'); l3 = 'Traditional - AVG RE';
%                         end
                    end
                    ylabel('Relative Error')
                    xlabel('Number of snapshots/function counts per injection pattern')
%                     if ~isempty(find(obj.idx,1))
%                         legend(l1,l2,l3,l4);
%                     else
%                         legend(l1);
%                     end
                    legend('FontSize',15,'Location','best');
                    grid();
                    set(gca,'FontSize',15);
                end
                
            elseif ~isempty(obj.split_conds)

                for i = 1:obj.num_samples
                    semilogy(1:obj.range,obj.errors_ROM_split(:,:,i),'-^','color','r'); l1 = 'ROM';
                    hold on
                    if ~isempty(find(obj.idx,1))
                        semilogy(1:obj.range,obj.errors_TRAD_split(:,:,i),'-o','color','k'); l2 = 'TRAD';
                    end
                end
                ylabel('Relative Error in conductivity estimation')
                xlabel('Number of snapshots/function counts per injection pattern')
                if ~isempty(find(obj.idx,1))
                    legend(l1,l2);
                else
                    legend(l1);
                end
                grid();
                
            end
        end
        
        function obj = saveSnap(obj)
            disp('Saving processed Snapshot comparison data in ROM/Results/snap_data.mat')
            snap=obj;
            OrderedModelClass.changePath('ROM'); top = getenv("ROMEG_TOP");
            save([top '/Results/snap_data.mat'],'snap')
        end
    end 
end