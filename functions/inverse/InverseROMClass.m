classdef InverseROMClass < InverseClass & OrderedModelClass

    properties
        LF              % ROM model
        simultaneous    % simultaneous electrode estimation
        snaps           % boolean, run inverse for all number of snapshots
        snap            % current number of snapshots to use in RBModel
        min_snap
        tag = ''        % how to label the estimates
        weights
        weighted
        omit_layers
    end

    methods

        function obj = InverseROMClass(varargin)

            %obj = obj.processArgs(varargin);
            obj@OrderedModelClass(varargin)
            obj = obj.getTOP();
            obj = obj.loadSinks();
            if ~isempty(obj.weighted) && obj.weighted
                obj = obj.loadWeights();
                obj.logger.info('InverseROMClass','Loading measurement weightings')
            end

        end

        function obj = setUp(obj)
            
            obj.el_in = obj.sinks(obj.pattern,1);
            obj = obj.loadLF();
            
            LF = obj.LF{obj.el_in};
            
            if isempty(obj.active_layers)
                obj.lf = LF.non_active;
                obj.te = LF.active;
                obj.active_layers = LF.active;
                disp('Active layers set as the number used for training.')
                disp('To change this use the active_layers option in the GenInverse function.')
            else
                non_active_layers = 1:LF.P;
                non_active_layers(obj.active_layers) = [];
                obj.lf = non_active_layers;
                obj.te = obj.active_layers;
            end

            if obj.fix_conds
                obj.cond_lf = (LF.mu_min(obj.lf) + LF.mu_max(obj.lf))/2;
            else
                if isempty(obj.omit_layers)
                    obj.cond_lf = obj.synth_cond(obj.lf);
                else
                    conds_tmp = obj.synth_cond;
                    conds_tmp(obj.omit_layers) = [];
                    obj.cond_lf = conds_tmp(obj.lf);
                end
                %disp('Warning: Using synth conds from measurements to assign layers')
                %disp('Make sure they are the same length and if not use fix_conds')
            end
            obj.lb = LF.mu_min(obj.te);
            obj.ub = LF.mu_max(obj.te);
            
            if isempty(obj.x0), obj.x0 = (obj.lb + obj.ub)/2; end

        end

        function obj = opt(obj)
            
            if ~isempty(obj.simultaneous) && obj.simultaneous
                options = optimoptions(@fmincon,'Display','iter','Algorithm','interior-point',...
                    'FiniteDifferenceType','central','OptimalityTolerance',1e-13,...
                    'MaxFunctionEvaluations',5000,'MaxIterations',2000);
            else
                options = optimoptions(@fmincon,'Display','iter','Algorithm','interior-point',...
                    'FiniteDifferenceType','central','OptimalityTolerance',1e-13,...
                    'MaxFunctionEvaluations',20000,'MaxIterations',2000);
            end
            
            if ~isempty(obj.snaps) && obj.snaps
                runs = 1:obj.min_snap;
            else
                runs = 1;
            end
            
            obj.estimate = [];
            for i = runs
                if ~isempty(obj.snaps) && obj.snaps
                    obj.snap = i;
                else
                    obj.snap = [];
                end
                
                func=@(cond_te)obj.functionEITSim(cond_te);
                estimate=fmincon(func,obj.x0,obj.A,obj.b,obj.Aeq,obj.beq,obj.lb,obj.ub,obj.nonlcon,options);
                obj.estimate = [obj.estimate; estimate];
                obj.logger.info('opt',['The estimated conductivities are ' num2str(estimate)])
                obj.logger.info('opt',['The synthetic conductivities are ' num2str(obj.synth_cond)])
            end
        end
        
        function zN = RBsolution(obj,num,mu_a)
            % returns the potential given by RB solution on the electrodes
            
            if isempty(obj.snap)
                obj.snap = size(obj.LF{num}.V,2);
            end            
            
            n_mu=length(mu_a); % number of parameters
            M_mu=obj.LF{num}.ANq{n_mu}(1:obj.snap,1:obj.snap)/mu_a(n_mu); % Z
            
            for kk=1:n_mu-1
                M_mu=M_mu+mu_a(kk)*obj.LF{num}.ANq{kk}(1:obj.snap,1:obj.snap); % Stiffness
            end
            
            zN=M_mu\obj.LF{num}.FNq{1}(1:obj.snap);
            
        end
       
	function [zNh,zN] = RBapprox(obj,num,mu_a)
	    
	    zN = obj.RBsolution(num,mu_a);

	    if isempty(obj.snap)
            obj.snap = size(obj.LF{num}.V,2);
	    end

	    zNh = obj.LF{num}.V(:,1:obj.snap)*zN;
	end

        function zNh = combinedRBsolution(obj,mu_a,el_in,el_out)            
            
            if ~isempty(obj.new_sinks) && obj.new_sinks %&& isempty(obj.use_sinks)
                [zNh,~] = obj.RBapprox(el_in,mu_a);
                for jj = el_out
                    if ~(jj == obj.ref_sink)
                        [z_tmp,~] = obj.RBapprox(jj,mu_a);
                        zNh = zNh - z_tmp/length(el_out);
                    end
                end
                zNh = zNh(end-(obj.eL-1):end);
                zNh([el_in el_out]) = [];
            else
                [zNh,~] = obj.RBapprox(el_in,mu_a);
                zNh = zNh(end-(obj.eL-1):end);
                zNh([el_in el_out]) = [];
            end
            
        end
        
        function mu_a = makeMu(obj,cond_te)
            
            N_layers=length([obj.te obj.lf]); % Number of layers to estimate
            lis=1:N_layers;lis(obj.lf)=[];
            
            mu_a = zeros(1,N_layers);
            for ii=1:length(obj.lf)
                mu_a(obj.lf(ii)) = obj.cond_lf(ii);
            end
            
            for ii=1:length(lis)
                mu_a(lis(ii)) = cond_te(ii);
            end
        end
        
        function f = functionEITSim(obj,cond_te)

            mu_a = obj.makeMu(cond_te);
            f_tmp = zeros(obj.eL,1);
            if ~isempty(obj.simultaneous) && obj.simultaneous
                for ii = 1:size(obj.sinks,1)
                    %obj.pattern = ii;
                    el_in = obj.sinks(ii,1); el_out = obj.sinks(ii,2:end);
                    zNh1 = obj.combinedRBsolution(mu_a,el_in,el_out);

                    % Compute error between measurement and simulation
                    %disp(size(obj.u{ii}))
                    u = obj.u{ii}(end -(obj.eL-1):end);
                    u([el_in el_out]) = [];
                    di=zNh1-u;
                    if ~isempty(obj.weighted) && obj.weighted
                        %weights = normalize(obj.weights(ii,:)',"norm",1);
                        %obj.weights([el_in el_out]) = [];
                        di = di.*obj.weights(ii,:)';
                    end
                    f_tmp(ii,1)=norm(di);%/norm(zNh1);
                end
                f = mean(f_tmp,1);
            else
                el_in = obj.sinks(obj.pattern,1); el_out = obj.sinks(obj.pattern,2:end);
                zNh1 = obj.combinedRBsolution(mu_a,el_in,el_out);
                
                u = obj.u{obj.pattern}(end -(obj.eL-1):end);
                u([el_in el_out]) = [];
                di=zNh1-u;
                if ~isempty(obj.weighted) && obj.weighted
                    %weights = normalize(obj.weights(obj.pattern,:)',"norm",1);
                    %obj.weights([el_in el_out]) = [];
                    di = di.*obj.weights(obj.pattern,:)';
                end
                f=norm(di)/norm(zNh1);
            end 
        end

        function obj = loadLF(obj)
            if isempty(obj.LF)
                obj.logger.info('loadLF','Loading RBModel from Results folder')
                load([obj.top '/Results/ROM/RBModel.mat'],'RBModel')
            end
            if ~isempty(obj.new_sinks) && obj.new_sinks && isempty(obj.simultaneous) %&& isempty(obj.use_sinks)
                obj.LF = {};
                N_list = [];
                sink_elec = RBModel.LF{1}.sink_elec;
                for i=obj.sinks(obj.pattern,:)
                    if ~(i==sink_elec)
                        %disp(['Loading LF for pattern ' num2str(i)])
                        obj.LF{i} = RBModel.LF{i};
                        obj.logger.debug('loadLF',['Loading ROM LF ' num2str(i) ' into obj.LF ' num2str(i)])
                        N_list = [N_list RBModel.LF{i}.N];
                    end
                end
                obj.eL = RBModel.LF{obj.el_in}.L;
                obj.min_snap = min(N_list);
            elseif ~isempty(obj.simultaneous) && (obj.simultaneous == true)
                obj.LF = RBModel.LF;
                obj.logger.debug('loadLF','Loading full ROM LF into obj.LF')
                obj.eL = RBModel.LF{1}.L;
                obj.min_snap = RBModel.LF{1}.N;
            else
                obj.LF{obj.el_in} = RBModel.LF{obj.pattern};
                obj.logger.debug('loadLF',['Loading ROM LF ' num2str(obj.pattern) ' into obj.LF ' num2str(obj.el_in)])
                obj.eL = RBModel.LF{obj.pattern}.L;
                obj.min_snap = RBModel.LF{obj.pattern}.N;
            end
        end

        function saveInv(obj)
            
            folder = 'inverse_';
            for ii = 1:length(obj.active_layers)
                folder = [folder num2str(obj.active_layers(ii))];
            end
            
            inv = struct();
            inv.estimate = obj.estimate;
            inv.synth_cond = obj.synth_cond;
            inv.cond_lf = obj.cond_lf;
            save([obj.top '/Results/inverse/ROM/' folder '/inv_' num2str(obj.pattern) '.mat'], 'inv')
        end

        function savePrep(obj)
            inv = obj;
            save([obj.top '/Results/inverse/ROM/prep.mat'],'inv')
        end

        function obj = collect(obj)
            
            folder = 'inverse_';
            for ii = 1:length(obj.active_layers)
                folder = [folder num2str(obj.active_layers(ii))];
            end
            
            if ~isempty(obj.simultaneous) && obj.simultaneous
                obj.num_patterns = 1;
            end
            
            if obj.snaps
                estimates = [];
                for i=1:obj.num_patterns
                    load([obj.top '/Results/inverse/ROM/' folder '/inv_' num2str(i) '.mat'],'inv')
                    estimates{i} = inv.estimate;
                    delete([obj.top '/Results/inverse/ROM/' folder '/inv_' num2str(i) '.mat'])
                end
                %estimate = mean(estimates,3);
                sinks = obj.sinks;
                save([obj.top '/Results/inverse/ROM/' folder '/' obj.tag '_estimate_snaps.mat'],'estimates','sinks')
                disp(['Collected electrode estimates and saved result (with sinks) in Results/inverse/ROM/' folder '/' obj.tag '_estimate_snaps.mat'])
            else
                estimates = [];
                for i=1:obj.num_patterns
                    load([obj.top '/Results/inverse/ROM/' folder '/inv_' num2str(i) '.mat'],'inv')
                    estimates = [estimates; inv.estimate];
                    delete([obj.top '/Results/inverse/ROM/' folder '/inv_' num2str(i) '.mat'])
                end
                obj.estimates = estimates;
                estimate = mean(estimates,1);
                sinks = obj.sinks;
                disp(folder)
                save([obj.top '/Results/inverse/ROM/' folder '/' obj.tag '_estimate.mat'],'estimate','sinks')
                save([obj.top '/Results/inverse/ROM/' folder '/' obj.tag '_estimates.mat'],'estimates','sinks')
                disp(['The average estimated conductivity using the ROM method is ' num2str(estimate)])
                disp(['The synth conductivity values ---------------------------> ' num2str(inv.synth_cond)])
                disp(['The fixed conductivity values ---------------------------> ' num2str(inv.cond_lf)])
                disp(['Collected electrode estimates, averaged them and saved result (with sinks) in Results/inverse/ROM/' folder '/' obj.tag '_estimate.mat'])
            end
        end
        
        function obj = loadWeights(obj)
          
            obj.logger.info('loadLF','Loading weights from ROM folder')
            load([obj.top '/Results/ROM/weights.mat'],'weights')
            obj.weights = weights;
        end
    end
end
