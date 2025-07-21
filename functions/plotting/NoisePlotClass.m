classdef NoisePlotClass < PlottingClass
    properties
        errors_ROM
        avg_RE
        noise_levels
        snapshot
        layer_names
    end
    
    methods
        
        function obj = NoisePlotClass(varargin)
        %
        %   NoisePlotClass(name1,value1,name2,value2,...)
        %
        % Arguments:
        %   range           - (int) visualize how many snapshots
        %
        %
            obj = obj.processArgs(varargin);
            
        
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
                for j = 1:length(obj.noise_levels)
                    tag2 = [obj.tag num2str(j)];
                    load([obj.top '/Results/inverse/ROM/' folder2 '/' tag2 '_estimate_snaps.mat'],'estimates')
                    for k=1:length(estimates)
                        if size(estimates{k},1) < 100
                            estimates{k} = [estimates{k}; NaN(100-size(estimates{k},1),size(estimates{k},2))];
                        end
                        estimates1(:,:,k) = estimates{k};
                    end
                    estimates2 = mean(estimates1,3,"omitnan");
                    obj.results_ROM{j}(:,:,i) = estimates2;
                end
                load([obj.top '/Results/measurements/prep.mat'],'Data')
                layers = c;
                obj.conds(:,:,i) = Data.synth_cond(layers);
            end
        end
        
        function obj = processROMResults(obj)
            for j = 1:length(obj.noise_levels)
                samples = size(obj.results_ROM{j},3);
                for i = 1:length(samples)
                    obj.errors_ROM{j} = abs(obj.results_ROM{j}(:,:,i)-obj.conds(:,:,i))./obj.conds(:,:,i);
                end
                obj.avg_RE(:,:,j) = mean(obj.errors_ROM{j},3,"omitnan");
            end
        end
        
        function plotNoise(obj,varargin)
            
            figure()
            
            icons = [{'-diamond'},{'-square'},{'-*'},{'-^'},{'-o'},{'-x'}];
            colors = [{"#0072BD"},{'r'},{'k'},{"#D95319"},{"#EDB120"},{"#7E2F8E"}];
            for i = obj.tissue
                data = reshape(obj.avg_RE(obj.snapshot,i,:),1,[]);
                semilogy(obj.noise_levels,data,icons{i},'color',colors{i},'DisplayName',obj.layer_names{i});
                hold on
            end
            ylabel('RE')
            xlabel('\textup{Standard Deviation in Gaussian Noise $(\mu V)$}','interpreter','latex','FontSize',18)
            %\fontname{Arial} Number of \it{n} \times \it{n}\rm linear systems to solve','interpreter','tex','FontSize',15

            legend('FontSize',15,'Location','best');
            grid();
            set(gca,'FontSize',15);
        end
    end 
end