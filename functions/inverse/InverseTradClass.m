classdef InverseTradClass < InverseClass

    properties
        S           % Prelim stiffness matrix
        source      % source vector
        history
        snaps
    end

    methods

        function obj = InverseTradClass(varargin)
            obj = obj.processArgs(varargin);
            obj = obj.loadSinks();
            %obj = obj.processModel();
        end

        function obj = setUp(obj)

            obj.el_in = obj.sinks(obj.injection,1);
            el_in = obj.el_in; el_out = obj.sinks(obj.injection,2:end);

            FOM = FOMClass.loadFOM(obj.top,false);
            p = FOM.p; t = FOM.t; f = FOM.f; 
            
            obj.Ind_E = zeros(FOM.L,1);
            for ii = 1:FOM.L
                indx = unique(f(f(:,4)==ii,1:3));
                obj.Ind_E(ii,1) = indx(round(length(indx)/2));
            end

            obj.eL = FOM.L; np = FOM.np;

%             obj.lf = FOM.non_active; obj.cond_lf = obj.synth_cond(FOM.non_active);
%             obj.lb = FOM.mu_min(FOM.active);
%             obj.ub = FOM.mu_max(FOM.active);
%             obj.te = FOM.active;
            
            if isempty(obj.active_layers)
                obj.lf = FOM.non_active;
                obj.te = FOM.active;
                obj.active_layers = FOM.active;
                disp('Active layers set as the number used for training.')
                disp('To change this use the active_layers option in the GenInverse function.')
            else
                non_active_layers = 1:FOM.P;
                non_active_layers(obj.active_layers) = [];
                obj.lf = non_active_layers;
                obj.te = obj.active_layers;
            end

            if obj.fix_conds
                obj.cond_lf = (FOM.mu_min(obj.lf) + FOM.mu_max(obj.lf))/2;
            else
                obj.cond_lf = obj.synth_cond(obj.lf);
            end
            obj.lb = FOM.mu_min(obj.te);
            obj.ub = FOM.mu_max(obj.te);
            
            if isempty(obj.x0), obj.x0 = (obj.lb + obj.ub)/2; end

            I_in = obj.current;
            I_out=(-1*I_in)/length(el_out);
            obj.source = spalloc(size(p,1)+obj.eL,1,length(el_out)+1);
            obj.source(np+el_in) = I_in;
            obj.source(el_out+np) = I_out;
            
            N_stiff = length([obj.lf obj.te]);

            obj.S=prepare_EIT_model_cem(p,t,obj.eL,np);
            [A,B,C] = femeg_stiffness_cem(p,t,f,1);
            S_cem = [A, -B; -B', C];
            obj.S(N_stiff).stiff = S_cem;
            
            disp(obj)

        end

        function obj = opt(obj)

            history = [];
            func=@(cond_te)obj.functionEIT(cond_te);
            options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton',...
                'SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-13,'TolX',1e-12,...
                'OutputFcn', @outputFunc);
            obj.estimate=fminunc(func,obj.x0,options);
            disp(['The estimated conductivity is ' num2str(obj.estimate)])

            obj.history = history; %#ok<*PROP> 

            function stop = outputFunc(x,optimValues,state)
    
            stop = false;
            
            switch state
                case 'init'
                    history = struct();
                    history.x = x;
                    history.funccount = optimValues.funccount;
                    history.fval = optimValues.fval;
                    history.iteration = optimValues.iteration;
                case 'iter'
                    ind = length(history)+1;
                    history(ind).x = x;
                    history(ind).funccount = optimValues.funccount;
                    history(ind).fval = optimValues.fval;
                    history(ind).iteration = optimValues.iteration;
                case 'done'
                otherwise
            end
            end
        end

        function [f,g] = functionEIT(obj,cond_te)
            tic
            N_layers=length([obj.te obj.lf]); % Number of layers to estimate
            lis=1:N_layers;lis(obj.lf)=[];
            
            mu_a = zeros(1,N_layers);
            for ii=1:length(obj.lf)
                mu_a(obj.lf(ii)) = obj.cond_lf(ii);
            end
            
            for ii=1:length(lis)
                mu_a(lis(ii)) = cond_te(ii);
            end
            
            % Compute stiffness matrix for a given conductivity
            n_mu=length(mu_a); % number of parameters
            Sn=obj.S(n_mu).stiff/mu_a(n_mu); % Z
            
            for kk=1:n_mu-1
                Sn=Sn+mu_a(kk)*obj.S(kk).stiff; % Stiffness
            end
            
            % Compute preconditioners
            [L,U]=ilu(Sn);
            [L1,U1] = ilu(Sn(1:end-obj.eL,1:end-obj.eL));
            
            % Solve FWD-P
            [un,~] = pcg(Sn,obj.source,1e-10,6000,L,U);
            un2 = un(end-(obj.eL-1):end);
            disp(['The time after full order calc is ' num2str(toc)])
            %un2 = un(Ind_E); un2(el_in) = []; un2(el_out) = [];
            %u22 = u2(Ind_E); u22(el_in) = []; u22(el_out) = [];
            
            % Compute error between measurement and simulation
            di=un2-obj.u(end-(obj.eL-1):end);
            %di = di(Ind_E);
            f=norm(di);
            
            % Compute gradient
            g=zeros(length(obj.Ind_E),length(lis));cont=1;
            for kk=lis
                [gg,~] = pcg(Sn(1:end-obj.eL,1:end-obj.eL),obj.S(kk).stiff(1:end-obj.eL,1:end-obj.eL)*un(1:end-obj.eL),1e-10,6000,L1,U1);
                g(:,cont)=gg(obj.Ind_E); %MK^-1K_iv
                cont=cont+1;
            end
            %g([el_in,el_out],:)=[];
            g=-2*di'*g;
            disp(['The time after gradient calc is ' num2str(toc)])
        end

        function saveInv(obj)
            
            folder = 'inverse_';
            for ii = 1:length(obj.active_layers)
                folder = [folder num2str(obj.active_layers(ii))];
            end
            
            inv = obj;
            inv.S = [];
            inv.source = [];
            inv.u = [];
            save([obj.top '/Results/inverse/TRAD/' folder '/inv_' num2str(obj.injection) '.mat'], 'inv')
        end

        function savePrep(obj)
            inv = obj;
            save([obj.top '/Results/inverse/TRAD/prep.mat'],'inv')
        end
        
        function collect(obj)
            
            folder = 'inverse_';
            for ii = 1:length(obj.active_layers)
                folder = [folder num2str(obj.active_layers(ii))];
            end
            
            estimates = [];
            histories = [];
            for i=1:obj.num_patterns
                if isfile([obj.top '/Results/inverse/TRAD/' folder '/inv_' num2str(i) '.mat'])
                    load([obj.top '/Results/inverse/TRAD/' folder '/inv_' num2str(i) '.mat'],'inv')
                    estimates = [estimates; inv.estimate];
                    histories{i} = inv.history;
                    delete([obj.top '/Results/inverse/TRAD/' folder '/inv_' num2str(i) '.mat'])
                else
                    estimate_empty = zeros(1,length(obj.active_layers));
                    estimates = [estimates; estimate_empty];
                    histories{i}=[];
                end
            end
            estimate = mean(estimates,1);
            save([obj.top '/Results/inverse/TRAD/' folder '/estimate.mat'],'estimates','histories')
            disp(['The average estimated conductivity using the traditional method is ' num2str(estimate)])
            disp('Collected electrode estimates, averaged them and saved them in Results/inverse/TRAD/estimate.mat')
        end
    end
end