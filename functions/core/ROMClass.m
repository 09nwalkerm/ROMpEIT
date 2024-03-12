classdef ROMClass < OrderedModelClass
    properties
        LF
        N
        V
        Qa
        Qf
        P
        L
        active
        non_active
        mu_max
        mu_min
        sample_grid
        Greedy_samples
        ANq
        FNq
        Xnorm
        delta_Max
        delta_Mean
        delta
        tolGREEDY = 5e-5    % default value for tolGREEDY = 5e-5
        Nmax = 30           % default value for Nmax=30
        Cqq
        dqq
        Eqq
        XAV
        Z
        FOM
        electrode
        injection
        sink_elec
        current             % injection current
        L1
        U1
        tol
        split
    end
    methods
        function obj = ROMClass(varargin)
        %
        % ROMClass(name1, value1, name2, value2 ...)
        %
        % Description:
        %   A class object containing the information and functions needed
        %   to produce the RBModel. Please use the GenRBModel function to 
        %   generate the RBModel. Although feel free to use this constructor
        %   function for debugging purposes or if using a pre-calculated 
        %   stability factor interpolant. 
        %
        % Args:
        %   electrode: (essential) The number of the desired injection
        %       electrode
        %   current: The injection current in Amps.
        %   tolGREEDY: The minumum error tolerance the greedy algorithm
        %       will stop at.
        %   Nmax: The maximum number of snapshots. If specified
        %       with tolGREEDY, the algorithm will stop at which 
        %       ever value is reached first.
        %   top: Path to top of ROMEG tree
        %
            
            %obj = obj.processArgs(varargin);
            obj@OrderedModelClass(varargin)
            obj = obj.getTOP();
            obj.FOM = FOMClass.loadFOM(obj.top);
            obj = obj.processArgs(obj.FOM.paramsROM);
            obj = obj.startLogger();
        end

        function obj = buildROM(obj)
            
            obj = obj.popFields();

            [obj,obj.FOM] = obj.solvePrelim(obj.FOM);

            obj = obj.sampleSpace(obj.FOM);

            obj = obj.runGreedy(obj.FOM);

            obj = obj.calcXnorm(obj.FOM);

            obj = obj.cleanROM();

        end

        function obj = calcXnorm(obj,FOM)
            obj.Xnorm = obj.V'*(FOM.Xnorm*obj.V);
        end

        function obj = popFields(obj)
            % Populate the ROMClass object fields from the arguments.
            
%             obj.FOM = FOMClass.loadFOM(obj.top);
%             
%             obj = obj.processArgs(FOM.paramsROM);

            if isempty(obj.current)
                obj.current = 0.020e-3;
            end
            
            if obj.use_sinks
                if ~isfield(obj,'top'), obj.top = getenv("ROMEG_TOP"); end
                obj = obj.loadSinks();
                obj.logger.info('popFields','Using sinks from file for injection and sink electrodes')
                obj.injection = obj.sinks(obj.electrode,1);
                obj.sink_elec = obj.sinks(obj.electrode,2:end);
            else
                obj.injection = obj.electrode;
                obj.sink_elec = obj.FOM.L;
            end

            obj.N=0;
            obj.V=[];
            obj.Qa=obj.FOM.Qa;
            obj.Qf=obj.FOM.Qf;
            obj.P=obj.FOM.P;
            obj.active = obj.FOM.active;
            obj.non_active = obj.FOM.non_active;
            obj.mu_max = obj.FOM.mu_max;
            obj.mu_min = obj.FOM.mu_min;
            obj.L = obj.FOM.L;
            
            %%% Bin p,t,f to save memory
            obj.FOM.p=[]; obj.FOM.t=[]; obj.FOM.f=[];
        end

        function [obj,FOM] = solvePrelim(obj,FOM)
            %%%%% Electrode matrix
            FOM.Fq{1}=sparse([FOM.np+obj.injection FOM.np+obj.sink_elec],ones(1,size(obj.sink_elec,2)+1),[obj.current -ones(1,size(obj.sink_elec,2))*(obj.current)/size(obj.sink_elec,2)],FOM.np+FOM.L,1);
            %%%%% Solve preliminary systems (Cqq)
            obj.tol=1e-8;
            obj.logger.info('solvePrelim','\n ** Computing Cqq ** \n');
            [obj.L1,obj.U1]=ilu(FOM.Xnorm);
            for q1=1:FOM.Qf
                [t,flag_t,res_t]=pcg(FOM.Xnorm,FOM.Fq{q1},obj.tol,2000,obj.L1,obj.U1);
                obj.logger.info('solvePrelim',['   Done ' num2str(q1) ' out of ' num2str(FOM.Qf) '. Flag: ' num2str(flag_t) '. Res: ' num2str(res_t)])
                for q2=1:FOM.Qf
                    obj.Cqq{q1,q2}=t'*FOM.Fq{q2};
                end
            end
        end

        function obj = sampleSpace(obj,FOM)
            mu_train_Dimension=6000;
            mu_cube=lhsdesign(mu_train_Dimension,length(FOM.active)); % normalized design
            mu_train=bsxfun(@plus,FOM.mu_min(FOM.active),bsxfun(@times,mu_cube,(FOM.mu_max(FOM.active)-FOM.mu_min(FOM.active))));
            mu_1=(FOM.mu_max(FOM.active)+FOM.mu_min(FOM.active))/2;
            obj.sample_grid=[mu_1; mu_train];
        end

        function obj = runGreedy(obj,FOM)
            delta_Max=obj.tolGREEDY+1;
            new_mu_indx=1; M_sample = 1;
            sample_grid = obj.sample_grid;
            while obj.N < obj.Nmax && delta_Max > obj.tolGREEDY && M_sample > 1e-10
                tic
                    obj.N=obj.N + 1;
                    obj.logger.info('runGreedy',['\n **** Greedy iteration number ' num2str(obj.N) ' \n']);
                
                    mu_a=obj.sample_grid(new_mu_indx,:);
                    sample_grid(new_mu_indx,:)=[];
                    [M_mu,b_mu]=FOM.muAssemble(mu_a);[L_p,U_p]=ilu(M_mu);
                    [zh,flag_re]=pcg(M_mu,b_mu,1e-10,6000,L_p,U_p);disp(['Flag: ' num2str(flag_re)])

                    obj.logger.debug('runGreedy',['The time after full order calc is ' num2str(toc)])
                    zeta_n  = Gram_Schmidt_orth(obj.V, zh, FOM.Xnorm);
                    obj.logger.debug('runGreedy',['The time after gram-schmidt calc is ' num2str(toc)])
                    obj.V   = [obj.V zeta_n];
                    obj.Greedy_samples(obj.N,:) = mu_a;
                
                    [obj.ANq, obj.FNq] = project_System_EIT(FOM, obj.V);
                    obj.logger.debug('runGreedy',['The time after projection is ' num2str(toc)])
                    
                    [M_mu_N,~]=obj.muAssemble(ones(1,size(obj.sample_grid,2)));
                    M_sample = rcond(M_mu_N); % Used to test how singular the reduced stiffness matrix is becoming
                    obj.logger.debug('runGreedy',['Condition number of reduced stiffness matrix is currently ' num2str(M_sample)])
                    obj.logger.debug('runGreedy',['The time after rcond calc is ' num2str(toc)])
                    if (M_sample<1e-10)
                        obj.logger.info('runGreedy','Stopping Greedy Algorithm. Removing last snapshot...')
                        obj.V=obj.V(:,1:end-1);
                        obj.Greedy_samples = obj.Greedy_samples(1:end-1,:);
                        obj.N = obj.N -1;
                        [obj.ANq, obj.FNq] = project_System_EIT(FOM, obj.V);
                        obj.logger.info('runGreedy',['Reduced Order Model stopped at ' num2str(obj.N) ' snapshots.'])
                        toc
                        break
                    end
                
                    obj = femeg_offline_residual_iter_n(FOM,obj,obj.L1,obj.U1,obj.tol);
                    
                    obj.logger.debug('runGreedy',['The time after residual calc is ' num2str(toc)])
                    time1 = tic;
                    obj.logger.debug('runGreedy','Evaluate error estimate ... \n')
                    delta_N = zeros(1,size(sample_grid,1));

                    obj.Xnorm = obj.V(1:FOM.np,:)'*(FOM.Xnorm(1:FOM.np,1:FOM.np)*obj.V(1:FOM.np,:)); % re new
                    for ii = 1 : size(sample_grid,1)
                        
                        [M_mu_N,b_mu_N]=obj.muAssemble(sample_grid(ii,:));
                        zN=M_mu_N\b_mu_N;
                        
                        betaa=femeg_ROM_RBF_online(sample_grid(ii,:),FOM);betaa=real(betaa);

                        delta_N(ii)=femeg_ROM_error_estimate_eit(obj,zN,sample_grid(ii,:))/betaa/norm(zN);
                    end
                    
                    obj.logger.debug('runGreedy',['The time after error estimates is ' num2str(toc)])
                    obj.logger.debug('runGreedy',['Error estimates took: ' num2str(toc(time1))])
                    [delta_Max, new_mu_indx] = max(delta_N);
                    
                    delta_Mean=mean(delta_N);
                    obj.logger.info('runGreedy',[' Max Delta_N =  ' num2str(delta_Max)])
                    obj.logger.info('runGreedy',[' Mean Delta_N = ' num2str(delta_Mean)])
                
                    obj.delta_Max(obj.N) = delta_Max;
                    obj.delta_Mean(obj.N) = delta_Mean;
                    if length(obj.active)==1
                        obj.delta{obj.N} = delta_N;
                        obj.sample_grid{obj.N} = sample_grid';
                    end
%                     if obj.verbose
%                         ROM = obj;
%                         save([obj.top '/Results/verbose/snapshot_' num2str(obj.N) '.mat'],'ROM');
%                     end
                toc
            end
        end

        function obj = cleanROM(obj)
        % Some cleaning

            %obj.V=obj.V(end-obj.L+1:end,:);
            obj.dqq = [];
            obj.Eqq=[];
            obj.Cqq=[];
            obj.Z=[];
            obj.FOM=[];
            obj.L1 = [];
            obj.U1 = [];
            obj.XAV = [];
        end

        function saveLF(obj,electrode)
            ROM = obj;
            save([obj.top '/Results/ROM/other/LF_EIT_' num2str(electrode) '.mat'],'ROM',"-v7.3")
        end
    end
end
