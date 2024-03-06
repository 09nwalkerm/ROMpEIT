classdef SenseClass < PlottingClass
    properties
        ROM_std
        ROM_fullstd
        num_layers
        layer_names
        layers
        plot
        anis
        split_elec
        elec
        tag
        c
    end
    
    methods
        
        function obj = SenseClass(varargin)
        %
        %   SenseClass(name1,value1,name2,value2,...)
        %
        % Arguments:
        %   layer_names     - list of names in a cell
        %   layers          - layers to be estimated
        %   plot            - type of plot, "Bar", "Box"
        %   num_samples     - array of sample numbers to read from
        %   anis            - isotropic layer to be compared to
        %                     anisotropic counter parts
        %   split_elec      - split RE into individual electrodes
        %   c               - used to make folder string, layers estimated
        %   tag             - tag used to save estimates
        %
            
            obj = obj.processArgs(varargin);
            obj = obj.readResults();
            obj = obj.processResults();
            obj.plotSensitivity();
            
        end
        
        function obj = readResults(obj)

            if isempty(obj.num_samples)
                count=1; folder = ['Result' num2str(count)];
                OrderedModelClass.changePath(folder); obj.top = getenv("ROMEG_TOP");                
                while isfolder([obj.top '/../' folder])
                    OrderedModelClass.changePath(folder); obj.top = getenv("ROMEG_TOP");
                    for k = 1:size(obj.c,1)
                        folder2 = 'inverse_';
                        for l = 1:size(obj.c,2)
                            folder2 = [folder2 num2str(obj.c(k,l))];
                        end
                        load([obj.top '/Results/inverse/ROM/' folder2 '/' obj.tag '_estimate.mat'],'estimate')
                        obj.results_ROM(k,:,count) = estimate;
                    end
                    load([obj.top '/Results/measurements/prep.mat'],'Data')
                    %obj.layers=find(logical(Data.mu_max-Data.mu_min));
                    obj.conds(:,:,count) = Data.synth_cond;%(obj.layers);
                    obj.num_layers = length(Data.synth_cond);
                    count=count+1;
                    folder = ['Result' num2str(count)];
                end
            else
                for i = 1:length(obj.num_samples)
                    folder = ['Result' num2str(obj.num_samples(i))];
                    OrderedModelClass.changePath(folder); obj.top = getenv("ROMEG_TOP");
                    folder2 = 'inverse_';
                    for l = 1:size(obj.c,2)
                        folder2 = [folder2 num2str(obj.c(1,l))];
                    end
                    if isempty(obj.split_elec)
                        load([obj.top '/Results/inverse/ROM/' folder2 '/' obj.tag '_estimate.mat'],'estimate')
                        obj.results_ROM(:,:,i) = estimate; %TODO: split for patterns
                    else
                        load([obj.top '/Results/inverse/ROM/' folder2 '/' obj.tag '_estimates.mat'],'estimates')
                        obj.results_ROM(:,:,i) = estimates;
                    end
                    load([obj.top '/Results/measurements/prep.mat'],'Data')
                    %obj.layers=find(logical(Data.mu_max-Data.mu_min));
                    obj.conds(:,:,i) = Data.synth_cond;%(obj.layers);
                    obj.num_layers = size(obj.c,2);
                end
            end
        end
        
        function obj = processResults(obj)
            %obj.num_samples = size(obj.results_ROM,3);
            for ii = 1:length(obj.num_samples)
                for jj = 1:size(obj.c,1)
                    if ~isempty(obj.anis)
                        conds = [obj.conds(:,1:obj.anis,ii) obj.conds(:,obj.anis,ii) obj.conds(:,obj.anis+1:end,ii)];
                        conds = conds(:,obj.c(jj,:),:);
                    else
                        conds = obj.conds(:,obj.layers,ii);
                        %disp(conds)
                        conds = conds(:,obj.c(jj,:));
                        %conds = obj.conds(:,:,ii);
                    end                    

                    if isempty(obj.layers)
                        obj.ROM_RE(jj,:,ii) = abs(obj.results_ROM(jj,:,ii)-conds(obj.c(jj,:)))./conds(obj.c(jj,:));
                    else
                        %disp(jj);disp(ii);disp(conds)
                        obj.ROM_RE(jj,:,ii) = abs(obj.results_ROM(jj,:,ii)-conds)./conds;
                    end
                end
            end
            obj.ROM_std = mean(obj.ROM_RE,3);%std(obj.ROM_RE,0,3);
            obj.ROM_fullstd = NaN(size(obj.ROM_RE,1),obj.num_layers);
            for jj = 1:size(obj.ROM_RE,1)
                obj.ROM_fullstd(jj,obj.c) = obj.ROM_std(jj,:);
            end
        end
        
        function plotSensitivity(obj,varargin)
        %
        %   plotSensitivity(name1,value1,name2,value2,...)
        %
        % Arguments:
        %   layer_names     - list of names in a cell
        %   plot            - type of plot, "Bar", "Box"
        %   num_samples     - number of samples to read from
        %   anis            - isotropic layer to be compared to
        %                     anisotropic counter parts
        %   elec            - electrode to view sensitivity
        %
            obj = obj.processArgs(varargin);
            
            obj.num_layers = size(obj.layers,2);
            
            switch obj.plot
                case 'Bar'
                    obj.makeSensePlotBar(obj.layers)
                case 'Box'
                    obj.makeSensePlotBox(obj.layers)
                otherwise
                    obj.makeSensePlotBar(obj.layers)
            end
        end
        
        function makeSensePlotBar(obj,layers)
            
            figure()
            if isempty(obj.elec)
                y = obj.ROM_fullstd(1,:);
            else
                y = obj.ROM_fullstd(obj.elec,:);
            end
            x = categorical(obj.layer_names);
            x = reordercats(x,obj.layer_names);
            b = bar(x,y,'stacked','FaceColor','flat');
            ylabel('AVG RE of the estimation for single tissue')
            %b(1).DisplayName = 'Scalp';
            if (size(layers,1) > 1)
                for i = 1:size(layers,1)
                    for j = 1:size(layers,2)
                        cc=layers(i,:);
                        cc(:,j) = [];
                        z = cc;
                        disp(z)
                        b(i).CData(layers(i,j),:) = [0 1/z 1/z];
                    end
                end
            end
            hold on
            y = NaN(obj.num_layers);
            b2 = bar(x,y,'FaceColor','flat');
            if (size(layers,1) > 1)
                for i = 1:obj.num_layers
                    b2(i).DisplayName = obj.layer_names{i};
                    b2(i).FaceColor = [0 1/i 1/i];
                end
                legend(b2)
            end
            
        end
        
        function makeSensePlotBox(obj,layers)
            
            figure()
            if isempty(obj.elec)
                y = obj.ROM_RE(1,:,:);
            else
                y = obj.ROM_RE(obj.elec,:,:);
            end
            results_matrix = reshape(y,obj.num_layers,size(obj.num_samples,2),[])';
            boxplot(results_matrix,'Labels',obj.layer_names(layers));
            ylabel('RE of the estimation for single tissue','FontSize',20)
            %ylim([0,0.2])
            ax = gca;
            ax.YAxis.Scale ="log";
            grid()
        end
        
        function obj = saveSense(obj)
            disp('Saving processed Sensitivity analysis data in ROM/Results/sense_data.mat')
            sense=obj;
            OrderedModelClass.changePath('ROM'); top = getenv("ROMEG_TOP");
            save([top '/Results/sense_data.mat'],'sense')
        end
    end 
end