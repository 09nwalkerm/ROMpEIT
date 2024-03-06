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
        tolGREEDY
        Nmax
        Cqq
        dqq
        Eqq
        XAV
        Z
        FOM
        electrode
        injection
        sink_elec
        current
        L1
        U1
        tol
        verbose
        split
        use_sinks
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
        % Arguments:
        %   - electrode  - (essential) The number of the desired injection 
        %                  electrode
        %   - current    - The injection current in Amps.
        %   - FOM        - The full order model.
        %   - tolGREEDY  - The minumum error tolerance the greedy algorithm
        %                  will stop at.
        %   - Nmax       - The maximum number of snapshots. If specified
        %                  with tolGREEDY, the algorithm will stop at which 
        %                  ever value is reached first.
        %   - top        - Path to top of ROMEG tree
        %   - verbose    - true or false: if true returns additional values in
        %                  RBModel and prints debugging information.
        %
        %

            [obj,obj.FOM] = obj.popFields(varargin);
            
        end

        function obj = buildROM(obj)
            
            FOM = obj.FOM;

            [obj,FOM] = obj.solvePrelim(FOM);

            obj = obj.sampleSpace(FOM);

            obj = obj.runGreedy(FOM);

            obj = obj.calcXnorm(FOM);

            obj = obj.cleanROM();

        end

        function obj = calcXnorm(obj,FOM)
            if obj.verbose
                disp("The dimensions of ROM.V...")
                disp(size(obj.V))
                disp("The dimensions of FOM.Xnorm...")
                disp(size(FOM.Xnorm))
            end
            obj.Xnorm = obj.V'*(FOM.Xnorm*obj.V);
        end

        function [obj,FOM] = popFields(obj,args)
            % Populate the ROMClass object fields from the arguments.
            
            obj = obj.processArgs(args);
            
            if find(logical(strcmp('FOM',args)))
                FOM = obj.FOM;
                obj.FOM = [];
            else
                FOM = FOMClass.loadFOM(obj.top,obj.verbose);
            end
            
            obj = obj.processArgs(FOM.paramsROM);
            
            if isempty(obj.tolGREEDY)
                obj.tolGREEDY = 5e-5; % 0.5% Relative error
            end

            if isempty(obj.Nmax)
                if obj.verbose
                    obj.Nmax = 10; % Suitable Nmax for debugging
                else
                    obj.Nmax = 30; % Nmax example
                end
            end

            if isempty(obj.current)
                obj.current = 0.020e-3;
            end
            
            if obj.use_sinks
                if ~isfield(obj,'top'), obj.top = getenv("ROMEG_TOP"); end
                obj = obj.loadSinks();
                disp('Using sinks from file for injection and sink electrodes')
                obj.injection = obj.sinks(obj.electrode,1);
                obj.sink_elec = obj.sinks(obj.electrode,2:end);
            else
                obj.injection = obj.electrode;
                obj.sink_elec = FOM.L;
            end

            obj.N=0;
            obj.V=[];
            obj.Qa=FOM.Qa;
            obj.Qf=FOM.Qf;
            obj.P=FOM.P;
            obj.active = FOM.active;
            obj.non_active = FOM.non_active;
            obj.mu_max = FOM.mu_max;
            obj.mu_min = FOM.mu_min;
            obj.L = FOM.L;
            
            %%% Bin p,t,f to save memory
            FOM.p=[]; FOM.t=[]; FOM.f=[];
        end

        function [obj,FOM] = solvePrelim(obj,FOM)
            %%%%% Electrode matrix
            FOM.Fq{1}=sparse([FOM.np+obj.injection FOM.np+obj.sink_elec],ones(1,size(obj.sink_elec,2)+1),[obj.current -ones(1,size(obj.sink_elec,2))*(obj.current)/size(obj.sink_elec,2)],FOM.np+FOM.L,1);
            %%%%% Solve preliminary systems (Cqq)
            obj.tol=1e-8;
            fprintf('\n ** Computing Cqq ** \n');
            [obj.L1,obj.U1]=ilu(FOM.Xnorm);
            for q1=1:FOM.Qf
                [t,flag_t,res_t]=pcg(FOM.Xnorm,FOM.Fq{q1},obj.tol,2000,obj.L1,obj.U1);
                disp(['   Done ' num2str(q1) ' out of ' num2str(FOM.Qf) '. Flag: ' num2str(flag_t) '. Res: ' num2str(res_t)])
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
                    fprintf('\n **** Greedy iteration number %d \n', obj.N);
                
                    mu_a=obj.sample_grid(new_mu_indx,:);
                    sample_grid(new_mu_indx,:)=[];
                    [M_mu,b_mu]=FOM.muAssemble(mu_a);[L_p,U_p]=ilu(M_mu);
                    [zh,flag_re]=pcg(M_mu,b_mu,1e-10,6000,L_p,U_p);disp(['Flag: ' num2str(flag_re)])

                    disp(['The time after full order calc is ' num2str(toc)])
                    zeta_n  = Gram_Schmidt_orth(obj.V, zh, FOM.Xnorm);
                    disp(['The time after gram-schmidt calc is ' num2str(toc)])
                    obj.V   = [obj.V zeta_n];
                    obj.Greedy_samples(obj.N,:) = mu_a;
                
                    [obj.ANq, obj.FNq] = project_System_EIT(FOM, obj.V);
                    disp(['The time after projection is ' num2str(toc)])
                    
                    [M_mu_N,~]=obj.muAssemble(ones(1,size(obj.sample_grid,2)));
                    M_sample = rcond(M_mu_N); % Used to test how singular the reduced stiffness matrix is becoming
                    disp(['Condition number of reduced stiffness matrix is currently ' num2str(M_sample)])
                    disp(['The time after rcond calc is ' num2str(toc)])
                    if (M_sample<1e-10)
                        disp('Stopping Greedy Algorithm. Removing last snapshot...')
                        obj.V=obj.V(:,1:end-1);
                        obj.Greedy_samples = obj.Greedy_samples(1:end-1,:);
                        obj.N = obj.N -1;
                        [obj.ANq, obj.FNq] = project_System_EIT(FOM, obj.V);
                        disp(['Reduced Order Model stopped at ' num2str(obj.N) ' snapshots.'])
                        toc
                        break
                    end
                
                    obj = femeg_offline_residual_iter_n(FOM,obj,obj.L1,obj.U1,obj.tol);
                    
                    disp(['The time after residual calc is ' num2str(toc)])
                    time1 = tic;
                    fprintf('Evaluate error estimate ... \n')
                    delta_N = zeros(1,size(sample_grid,1));

                    obj.Xnorm = obj.V(1:FOM.np,:)'*(FOM.Xnorm(1:FOM.np,1:FOM.np)*obj.V(1:FOM.np,:)); % re new
                    for ii = 1 : size(sample_grid,1)
                        
                        [M_mu_N,b_mu_N]=obj.muAssemble(sample_grid(ii,:));
                        zN=M_mu_N\b_mu_N;
                        
                        betaa=femeg_ROM_RBF_online(sample_grid(ii,:),FOM);betaa=real(betaa);

                        delta_N(ii)=femeg_ROM_error_estimate_eit(obj,zN,sample_grid(ii,:))/betaa/norm(zN);
                    end
                    
                    disp(['The time after error estimates is ' num2str(toc)])
                    disp(['Error estimates took: ' num2str(toc(time1))])
                    [delta_Max, new_mu_indx] = max(delta_N);
                    
                    delta_Mean=mean(delta_N);
                    fprintf(' Max Delta_N = %2.3e,   ', delta_Max)
                    fprintf(' Mean Delta_N = %2.3e \n', delta_Mean)
                
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

            if obj.verbose
                disp("Keeping additonal info in RBModel.mat for debugging or analysis")
            else
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
        end

        function saveLF(obj,electrode)
            ROM = obj;
            save([obj.top '/Results/ROM/other/LF_EIT_' num2str(electrode) '.mat'],'ROM',"-v7.3")
        end
    end
end