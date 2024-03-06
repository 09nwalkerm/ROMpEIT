classdef OrderedModelClass
    properties
        p
        t
        f
        Ind_E           % Tetrahedron number related to electrode
        anis_rad
        anis_tan
        angles          % true if model contains theta values (default: false)
        model           % path to model
        theta
        elec_height
        top             % path to top of the ROMEG tree
        sinks           % injection and sink patterns for each electrode
        sinks_path      % path to the sink.mat file
        num_sinks       % number of current sinks in the injection pattern
        num_patterns    % number of injection patterns
        electrodes      % list of injection electrodes for sink patterns
        out             % extraction electrodes
        new_sinks       % (boolean) use a new set of sinks called new_sinks.mat
        sample_num      % sample number to use for the inverse problem
        RB_path
        logLevel = 3    % level of logging for the log files, defualt is info
        cmdLevel = 3    % level of output for cmd logging, default is info
        num_dipoles
    end
    
    properties (Access = protected)
        logger
    end

    methods
        function obj = OrderedModelClass()
            obj.logger = log4m.getLogger();
            obj.logger.trace('OrderedModelClass','Loading logger into Class as obj.logger')
        end
        
        function [M_mu,b_mu] = muAssemble(obj,mu)
            
            if isa(obj,"ROMClass")
                matrix_field = 'ANq';
                source_field = 'FNq';
            elseif isa(obj,"FOMClass")
                matrix_field = 'Aq';
                source_field = 'Fq';
            end

            active = obj.active;
            non_active = obj.non_active;
            
            n_mu = zeros(1,obj.P);
            n_mu(active) = mu;
            n_mu(non_active) = obj.mu_min(non_active);
            
            M_mu=obj.(matrix_field){end}/n_mu(end);% Z

            for kk=1:length(n_mu)-1
                M_mu=M_mu+n_mu(kk)*obj.(matrix_field){kk}; % Stiffness
            end
            
            if nargout==2
                b_mu=obj.(source_field){1}(:,1);
            end
        end

        function obj = processArgs(obj,args)
            
            if ~isempty(args)
                
                if isa(args{1},'cell')
                    args = args{1};
                end

                for i = 1:2:length(args) % work for a list of name-value pairs
                    if ischar(args{i}) % check if is character
                        obj.(args{i}) = args{i+1}; % override or add parameters to structure.
                    end
                end
            end
        end

        function obj = processModel(obj)

            disp('Loading Model...')

            if isempty(obj.angles) || (~obj.angles)%isempty(obj.anis_rad) || (~isempty(obj.anis_rad) && (~obj.angles)) %|| ~isfield(obj,'anis_rad')
                try
                    load(obj.model,'p','t','f')
                    obj.p = p; obj.t = t; obj.f = f; %obj.Ind_E = Ind_E;
                    disp('Head model loaded')
                catch
                    fprintf("Cannot load head model, please check path given and that the file contains p,t,f,Ind_E values")
                    error('Cannot load head model')
                end
            else
                try
                    load(obj.model,'p','t','f','theta')
                    obj.p = p; obj.t = t; obj.f = f; obj.theta = theta;
                    disp('Head model loaded with angles')
                catch ME
                    switch ME.identifier
                        case 'MATLAB:load:couldNotReadFile'
                        fprintf("Cannot load head model, please check path given \n\n")
                        case 'MATLAB:UndefinedFunction'
                            fprintf("At least one variable missing, please ensure the head model file contains p,t,f,theta variables.\n\n")
                            fprintf("Or specify the name-value pair 'angles' false to have it generated.")
                            error('Missing variables')
                        otherwise
                            rethrow(ME)
                    end
                end
            end
        end

        function obj = loadSinks(obj,varargin)
            if isempty(obj.sinks)
                disp('Loading sinks from file...')
                if isempty(obj.sinks) && isempty(obj.sinks_path)
                    if obj.new_sinks
                        try
                            disp(['Loading sinks from ' obj.top '/Results/ROM/new_sinks.mat'])
                            load([obj.top '/Results/ROM/new_sinks.mat'],'sinks')
                        catch
                            disp('No path or new_sinks.mat file provided and none available in Results folder.')
                            disp("Please make new_sinks.mat using OrderedModelClass.patterns and rerun function.")
                            disp("Alternatively make one manually with custom injection patterns")
                            disp('See `help OrderedModelClass.patterns` for guidance.')
                            error('No new_sinks.mat file found')
                        end
                    else
                        try
                            disp(['Loading sinks from ' obj.top '/Results/ROM/sinks.mat'])
                            load([obj.top '/Results/ROM/sinks.mat'],'sinks')
                        catch
                            disp('No path or sinks.mat file provided and none available in Results folder.')
                            disp("Please make sinks.mat using OrderedModelClass.patterns and rerun function.")
                            disp("Alternatively make one manually with custom injection patterns")
                            disp('See `help OrderedModelClass.patterns` for guidance.')
                            error('No sinks.mat file found')
                        end
                    end
                elseif isempty(obj.sinks) && ~isempty(obj.sinks_path)
                    try
                        load(obj.sinks_path,'sinks')
                    catch
                        error('Cannot find .mat file with sinks. Path given could be incorrect.')
                    end
                end
                obj.sinks = sinks;
                obj.num_sinks = size(sinks,2)-1;
                obj.num_patterns = size(sinks,1);
            else
                obj.num_sinks = size(obj.sinks,2)-1;
                obj.num_patterns = size(obj.sinks,1);
            end
        end

        function obj = checkPaths(obj,varargin)
        % Arguments:
        %   type        - ROM, measurement, inverse
        %   num         - sample number (essential if
        %               type=measurement,eeg or inverse)
        %   RB_path     - path to RB model if NOT in $ROMEG_DATA
        %   num_dipoles - number of dipole folders to set up if using EEG
        %

            params = struct();
            for i = 1:2:length(varargin) % work for a list of name-value pairs
                if ischar(varargin{i}) % check if is character
                    params.(varargin{i}) = varargin{i+1}; % override or add parameters to structure.
                end
            end

            data = getenv("ROMEG_DATA");
            if strcmp(data(end),'/'), data = data(1:end-1); setenv("ROMEG_DATA",data); disp ('Resetting ROMEG_DATA'); end
            tree = getenv("ROMEG");
            if isfolder(data)
                % create log file if not already there
                if ~isfolder([data '/logs'])
                    mkdir([data '/logs'])
                end
                % create file system for results
                if ~strcmp(tree,data) || ~strcmp(tree,[data '/']) || ~strcmp([tree '/'],data)
                    if strcmp(params.type,'ROM')
                        if ~isfolder([data '/ROM'])
                            disp('Making ROM folder.')
                            mkdir([data '/ROM'])
                            OrderedModelClass.changePath('ROM');
                            obj.top = [data '/ROM'];
                            OrderedModelClass.setupFiles('ROM',true);
                        else
                            %disp('Warning: ROM folder in data path already exists, overwriting.')
                            OrderedModelClass.changePath('ROM')
                            obj.top = [data '/ROM'];
                        end
                    elseif strcmp(params.type,'measurement') || strcmp(params.type,'inverse') || strcmp(params.type,'bound')
                        if ~isfolder([data '/Result' num2str(params.num)])
                            disp(['Making Result' num2str(params.num) ' folder.'])
                            mkdir([data '/Result' num2str(params.num)])
                            OrderedModelClass.changePath(['Result' num2str(params.num)])
                            obj.top = [data '/Result' num2str(params.num)];
                            if isfield(params,'RB_path')
                                OrderedModelClass.setupFiles('ROM',false,'RB_path',params.RB_path)
                            else
                                OrderedModelClass.setupFiles('ROM',false)
                            end
                        else
                            disp(['Warning: Result' num2str(params.num) ' folder in data path already exists, overwriting.'])
                            OrderedModelClass.changePath(['Result' num2str(params.num)])
                            obj.top = [data '/Result' num2str(params.num)];

                        end
                    elseif strcmp(params.type,'eeg')
                        if ~isfolder([data '/Result' num2str(params.num)])
                            disp(['Making Result' num2str(params.num) ' folder.'])
                            mkdir([data '/Result' num2str(params.num)])
                            OrderedModelClass.changePath(['Result' num2str(params.num)])
                            obj.top = [data '/Result' num2str(params.num)];
                            OrderedModelClass.EEGFiles('sample_num',params.num,'num_dipoles',params.num_dipoles);
                        else
                            disp(['Warning: Result' num2str(params.num) ' folder in data path already exists, overwriting.'])
                            OrderedModelClass.changePath(['Result' num2str(params.num)])
                            obj.top = [data '/Result' num2str(params.num)];
                        end
                    end
                else
                    setenv("ROMEG_TOP",data)
                    obj.top = data;
                end
            end
        end

    end

    methods (Static)

        function patterns(varargin)
        %
        %   OrderedModelClass.patterns(name1,value1,name2,value2...)
        %
        %   Description:
        %       Function for making injection patterns with many sinks on
        %       the other side of the head model. Default behaviour is for
        %       every electrode but the last to be an injection and the
        %       last electrode being the extraction. If out is specified
        %       then all injection electrodes will have the extraction
        %       electrodes that are specified with the injection electrodes
        %       matching and of the extraction electordes being removed. 
        %       For opposite extraction electrodes use the num_sinks option.
        % 
        %   Arguments:
        %       model       - (essential) path to the head model
        %       num_sinks   - how many current sinks for each
        %                     injection would you like.
        %       elec_height - (optional) a given minimum height the sink
        %                     electrodes should be on the head model.
        %       electrodes  - list of electrodes for injection (defaults to
        %                     all)
        %       top         - (optional if ROMEG_TOP env var has been set)
        %                     path to the top of the ROMEG tree
        %       out         - specific list of extraction electrodes for
        %                     injection patterns.
        %       new_sinks   - are these new sinks being made for inverse
        %
        %   Examples:
        %       
        %       This example takes a model, with the electrode positions
        %       indicated in variable f and specifies 10 sink electrodes
        %       for each of the elctrodes listed, resulting in 8 patterns.
        %
        %       OrderedModelClass.patterns('model',model,'num_sinks',10,...
        %           'electrodes',[1,20,35,56,42,120,13,67])
        %

            obj = OrderedModelClass();
            obj = obj.processArgs(varargin);
            obj = obj.checkPaths('type','ROM');

            if isempty(obj.top)
                try
                    top = getenv("ROMEG_TOP");
                    obj.top=top;
                catch
                    disp('Please provide path to top of ROMEG tree using either "top" argument or by setting ROMEG_TOP environment variable')
                    error('Top of ROMEG tree is not defined')
                end
            end

            obj = obj.processModel();

            f = obj.f; p = obj.p; 

            disp('Generating sink patterns...')
            
            if isempty(obj.electrodes) && ~isempty(obj.num_sinks)
                obj.electrodes = 1:length(unique(f(:,4)))-1;
            elseif isempty(obj.electrodes) && isempty(obj.num_sinks)
                obj.electrodes = 1:length(unique(f(:,4)))-2;
            end

            if isempty(obj.out) && ~isempty(obj.num_sinks)
                elec_centers = zeros(length(unique(f(:,4)))-1,3);

                for ii = 1:length(elec_centers)
                    nodes = f(f(:,4)==ii,1:3);
                    nodes = unique(nodes);
                    elec_centers(ii,:) = [mean(p(nodes,1)) mean(p(nodes,2)) mean(p(nodes,3))];
                end

                sinks = zeros(length(obj.electrodes),obj.num_sinks+1);
                sinks(:,1) = obj.electrodes';

                for ii = 1:length(obj.electrodes)
                    pos = elec_centers(obj.electrodes(ii),:);

                    elec_sinks = elec_centers;
                    lengths = zeros(size(elec_sinks,1),3);

                    for jj = 1:size(elec_sinks,1)-1
                        lengths(jj,:) = elec_sinks(jj,:) - pos(:,:);
                        if ~isempty(obj.elec_height)
                            if (elec_sinks(jj,3) < obj.elec_height), lengths(jj,:) = []; end
                        end
                    end

                    norms = vecnorm(lengths');

                    [~,ind] = maxk(norms,obj.num_sinks);

                    sinks(ii,2:obj.num_sinks+1) = ind;
                end
            elseif ~isempty(obj.out)
                sinks = zeros(length(obj.electrodes),length(obj.out)+1);
                sinks(:,1) = obj.electrodes';
                sinks(:,2:end) = repmat(obj.out,length(obj.electrodes),1);
                sinks(obj.out,:) = [];
            else
                sinks = zeros(length(obj.electrodes),2);
                sinks(:,1) = obj.electrodes';
                sinks(:,2) = length(unique(f(:,4)))-1;
            end
            
            if obj.new_sinks
                save([obj.top '/Results/ROM/new_sinks.mat'],'sinks')
            else
                save([obj.top '/Results/ROM/sinks.mat'],'sinks')
            end
                
            disp('Saved sink patterns to /Results/ROM folder')
        end
        
        function wait(NAME,PAUSE,varargin)
            pause(5)
            if nargin > 2
                JOB=varargin{1};
            else
                [~,cmdout] = system(['squeue -u $USER | grep -m 1 ' NAME ' | awk ''{split($1,a,"_"); print a[1]; exit}'' ']);
                
                %[~,cmdout] = system(['squeue -u $USER | tail -n +3 | awk ''{split($1,a,"_"); print a[1]; exit}'' ']);
                cmdout2 = cmdout(1:end-1);
                JOB = cmdout2;
            end

%             while 1
%                 [~,cmdout] = system(['squeue -u $USER | grep ' NAME ' | awk ''{split($1,a,"_"); print a[1]; exit}'' ']);
%                 if isempty(cmdout)
%                     cmdout = cmdout(1:end-1);
%                     [~,cmdout2] = system(['scontrol show job ' cmdout ' | grep -B 3 ''JobState=FAILED'' | grep -w JobId | awk ''{split($1,a,"="); print a[2]; exit}'' ']);
%                     if isempty(cmdout2)
%                         break
%                     else
%                         while 1
%                             [~,cmdout3] = system(['scontrol show job ' cmdout ' | grep -B 3 ''JobState=FAILED'' | grep -w ArrayTaskId | awk ''{split($3,a,"="); print a[2]; exit}'' ']);
%                             if isempty(cmdout3)
%                                 break
%                             else
%                                 cmdout3 = cmdout3(1:end-1);
%                                 [CODE,~] = system(['scontrol requeue ' cmdout '_' cmdout3]);
%                                 if (CODE == 0); disp(['Requeued Job ' cmdout '_' cmdout3]); end
%                             end
%                         end
%                     end
%                 end
%                 pause(PAUSE)
%             end
            fprintf('\nWaiting for slurm job to finish')
            while 1
                fprintf('.')
                [~,cmdout] = system(['squeue --job ' JOB ' | grep -m 1 ' NAME ]);
                %[~,cmdout] = system(['squeue --job ' JOB ' | tail -n +2']);
                path = getenv("ROMEG");
                cmdofile = fopen([path '/scripts/cmdofile.txt'],'a+');
                fprintf(cmdofile,cmdout);
%                 [~,cmdout_tmp] = system(['which tail']);
%                 fprintf(cmdofile,cmdout_tmp);
%                 [~,cmdout_tmp] = system(['echo $PATH']);
%                 fprintf(cmdofile,cmdout_tmp);
                fprintf(cmdofile,cmdout2);
                fclose(cmdofile);
                if isempty(cmdout)
                    [~,cmdout2] = system(['scontrol show job ' JOB ' |grep -B 3 ''JobState=FAILED'' ']);
                    if isempty(cmdout2)
                        break
                    else
%                         while 1
%                             [~,cmdout3] = system(['scontrol show job ' JOB ' |grep -B 3 ''JobState=FAILED'' |grep -w ArrayTaskId |awk ''{split($3,a,"="); print a[2]; exit}'' ']);
%                             if isempty(cmdout3)
%                                 break
%                             else
%                                 cmdout3 = cmdout3(1:end-1);
%                                 [CODE,~] = system(['scontrol requeue ' JOB '_' cmdout3]);
%                                 if (CODE == 0); disp(['Requeued Job ' JOB '_' cmdout3]); end
%                             end
%                         end
                    end
                elseif strcmp(cmdout(1:end-1),'slurm_load_jobs error: Invalid job id specified')
                    disp(cmdout)
                    break
%                     [~,cmdout3] = system(['scontrol show job ' JOB ' |grep -B 3 ''JobState=FAILED'' |grep -w ArrayTaskId |awk ''{split($3,a,"="); print a[2]; exit}'' ']);
%                     if ~isempty(cmdout3)
%                         cmdout3 = cmdout3(1:end-1);
%                         [CODE,~] = system(['scontrol requeue ' JOB '_' cmdout3]);
%                         if (CODE == 0); disp(['Requeued Job ' JOB '_' cmdout3]); end
%                     end
                end
                pause(PAUSE)
            end
            
            while 1
                [~,cmdout] = system(['squeue -u $USER |grep ' NAME]);
                if isempty(cmdout),break;end
                [~,cmdout] = system(['squeue -u $USER |grep ' NAME ' |grep RH']);
                if ~isempty(cmdout)
                    [~,cmdout1] = system("squeue -u $USER |grep RH |awk '{print $1; exit}'");
                    [~,cmdout2] = system("squeue -u $USER |grep RH |awk '{split($1,a,""_""); print a[2]; exit}'");
                    [~,cmdout3] = system(['scancel ' cmdout1]);
                    setenv('JJOB',cmdout2)
                    setenv('NJOB',NAME);
                    !sbatch --array $JJOB -o $ROMEG_TOP/Results/slurm_logs/%x_%a_%j.out --job-name $NJOB $ROMEG/Functions/Cluster/cluster_job.sh $NJOB'
                end
                [~,cmdout] = system(['squeue -u $USER |grep ' NAME ' |grep PD |grep launch']);
                if ~isempty(cmdout)
                    [~,cmdout1] = system("squeue -u $USER |grep PD |grep launch |awk '{print $1; exit}'");
                    [~,cmdout3] = system(['scontrol hold ' cmdout1]);
                    pause(3)
                    [~,cmdout4] = system(['scontrol release ' cmdout1]);
                end
                pause(PAUSE)
            end
            fprintf('\n\n')
            %delete([path '/scripts/cmdofile.txt'])
        end

        function setupFiles(varargin)
        %
        %   setupFiles(name1,value1,name2,value2...)
        %
        % Description:
        %   Sets up folder structure for results from measurements and
        %   inverse runs. If not core then a symbolic link will be formed
        %   between the ROM folders.
        %
        % Arguments:
        %   ROM      - (boolean) will this folder contain the RBModel?
        %   RB_path     - path to RBModel for sym link if NOT in $ROMEG_DATA
        %
        %

            params = struct();
            for i = 1:2:length(varargin) % work for a list of name-value pairs
                if ischar(varargin{i}) % check if is character
                    params.(varargin{i}) = varargin{i+1}; % override or add parameters to structure.
                end
            end

            top = getenv("ROMEG_TOP");
            if ~isfolder([top '/Results'])
                disp('Setting up new Results folder')
                mkdir([top '/Results'])
                mkdir([top '/Results/verbose'])
                mkdir([top '/Results/slurm_logs'])
                mkdir([top '/Results/logs'])
                if ~params.ROM
                    if isfield(params,'RB_path')
                        setenv("ROMEG_RB",params.RB_path)
                        !ln -s $ROMEG_RB/ROM/Results/ROM %ROMEG_TOP/Results/ROM
                    else
                        !ln -s $ROMEG_DATA/ROM/Results/ROM $ROMEG_TOP/Results/ROM
                    end
                    mkdir([top '/Results/inverse'])
                    mkdir([top '/Results/bound'])
                    mkdir([top '/Results/inverse/ROM'])
                    mkdir([top '/Results/inverse/TRAD'])
                    mkdir([top '/Results/measurements'])
                else
                    mkdir([top '/Results/ROM'])
                    mkdir([top '/Results/ROM/other'])
                end
            else
                disp('Results folder already exists')
            end
        end

        function changePath(folder)
            %[~,path] = system(['readlink -f $ROMEG_DATA/' folder]);
            data = getenv("ROMEG_DATA");
            if strcmp(data(end),'/'), data = data(1:end-1); setenv("ROMEG_DATA",data); end
            path = [data '/' folder];
            if isfolder(path)
                setenv("ROMEG_TOP",path)
                disp(['ROMEG_TOP: ' path])
            else
                error('Folder not found')
            end
        end
        
        function sensitivityFiles(varargin)
        %
        %   OrderedModelClass.sensitivityFiles(name1,value1,name2,value2,...)
        %
        % Arguments:
        %   num_layers  - total number of active layers
        %   recursion   - number of layers to use
        %   sample_num  - sample number
        %   layers      - array of active layers where each row is new
        %                 active layer set
        %   order       - ROM or TRAD?
        %
        %
            
            params = struct();
            for i = 1:2:length(varargin) % work for a list of name-value pairs
                if ischar(varargin{i}) % check if is character
                    params.(varargin{i}) = varargin{i+1}; % override or add parameters to structure.
                end
            end
            
            if isempty(params.layers)
                layers = nchoosek(1:params.num_layers,params.recursion);
            else
                layers = params.layers;
            end
            
            OrderedModelClass.changePath(['Result' num2str(params.sample_num)])
            top = getenv("ROMEG_TOP");
            
            for i=1:size(layers,1)
                folder = 'inverse_';
                for j=1:size(layers,2)
                    folder = [folder num2str(layers(i,j))];
                end
                if ~isfolder([top '/Results/inverse/' params.order '/' folder])
                    mkdir([top '/Results/inverse/' params.order '/' folder])
                end
            end
        end
        
        function EEGFiles(varargin)
        %
        %   OrderedModelClass.EEGFiles(name1,value1,name2,value2,...)
        %
        % Arguments:
        %   sample_num  - sample number
        %   num_dipoles - number of dipoles being used
        %
            
            params = struct();
            for i = 1:2:length(varargin) % work for a list of name-value pairs
                if ischar(varargin{i}) % check if is character
                    params.(varargin{i}) = varargin{i+1}; % override or add parameters to structure.
                end
            end
            
            OrderedModelClass.changePath(['Result' num2str(params.sample_num)])
            top = getenv("ROMEG_TOP");
            if ~isfolder([top '/Results/EEG_FP'])
                mkdir([top '/Results/EEG_FP'])
                mkdir([top '/Results/slurm_logs'])
                mkdir([top '/Results/logs'])
            else
                warning('EEG_FP folder already exists, will be overwritten')
            end
            
            %folder = 'DP_';
            for i=1:params.num_dipoles
                folder = ['DP_' num2str(i)];
                if ~isfolder([top '/Results/EEG_FP/' folder])
                    mkdir([top '/Results/EEG_FP/' folder])
                end
            end
        end
        
        function backupROM(varargin)
        %
        %   OrderedModelClass.backupROM(name1,value1,name2,value2,...)
        %
        % Arguments:
        %   backup_folder  - name of the backup folder
        %   
        %
            params = struct();
            for i = 1:2:length(varargin) % work for a list of name-value pairs
                if ischar(varargin{i}) % check if is character
                    params.(varargin{i}) = varargin{i+1}; % override or add parameters to structure.
                end
            end
            
            data = getenv("ROMEG_DATA");
            if strcmp(data(end),'/'), data = data(1:end-1); setenv("ROMEG_DATA",data); end
            if ~isempty(params.backup_folder)
                if isfolder([data '/../' params.backup_folder])
                    copyfile([data '/ROM'],[data '/../' params.backup_folder '/ROM'])
                else
                    mkdir([data '/../' params.backup_folder])
                    copyfile([data '/ROM'],[data '/../' params.backup_folder '/ROM'])
                end
            end           
        end
        
        function resetInverse(varargin)
        %
        %   OrderedModelClass.resetInverse(name1,value1,name2,value2,...)
        %
        % Arguments:
        %   backup_folder   - name of the backup folder
        %   num_samples     - the number of samples to reset
        %
        %
        %
            params = struct();
            for i = 1:2:length(varargin) % work for a list of name-value pairs
                if ischar(varargin{i}) % check if is character
                    params.(varargin{i}) = varargin{i+1}; % override or add parameters to structure.
                end
            end

            data = getenv("ROMEG_DATA");
            if strcmp(data(end),'/'), data = data(1:end-1); setenv("ROMEG_DATA",data); end
            if ~isempty(params.backup_folder)
                if isfolder([data '/../' params.backup_folder])
                    for ii = 1:params.num_samples
                        OrderedModelClass.changePath(['Result' num2str(ii)])
                        top = getenv("ROMEG_TOP");
                        copyfile([top '/Results/inverse/ROM/inverse_*'],[data '/../' params.backup_folder '/Result' num2str(ii) '/Results/inverse/ROM/'])
                    end
                else
                    mkdir([data '/../' params.backup_folder])
                    for ii = 1:params.num_samples
                        OrderedModelClass.changePath(['Result' num2str(ii)])
                        top = getenv("ROMEG_TOP");
                        mkdir([data '/../' params.backup_folder '/Result' num2str(ii) '/Results/inverse/ROM/'])
                        copyfile([top '/Results/inverse/ROM/inverse_*'],[data '/../' params.backup_folder '/Result' num2str(ii) '/Results/inverse/ROM'])
                    end
                end
            end
            
            for ii = 1:params.num_samples
                OrderedModelClass.changePath(['Result' num2str(ii)])
                top = getenv("ROMEG_TOP");
                rmdir([top '/Results/inverse/ROM/inverse_*'] ,'s')
            end
        end
        
        function adjustROM()
            
            
            
            %load model
            OrderedModelClass.changePath('ROM')
            top = getenv("ROMEG_TOP");
            load([top '/Results/ROM/RBModel.mat'],'RBModel')
            
            %search and cut through V
            for ii=1:size(RBModel.LF,2)
                disp(['Cutting transformation matrix for pattern ' num2str(ii) ' of ' num2str(size(RBModel.LF,2))])
                condition = 1; count =1;
                while condition > 1e-10     % jj = 1:RBModel.LF{ii}.N
                    
                    m_mu = RBModel.LF{ii}.ANq{end}(1:count,1:count);
                    
                    for kk=1:RBModel.LF{ii}.P-1
                        m_mu = m_mu + RBModel.LF{ii}.ANq{kk}(1:count,1:count);
                    end
                    
                    condition = rcond(m_mu);
                    
                    if count == RBModel.LF{ii}.N
                        break
                    else
                        if condition < 1e-10
                            count = count - 1;
                            break
                        else
                            count = count + 1;
                        end
                    end
                end
                
                RBModel.LF{ii}.V = RBModel.LF{ii}.V(:,1:count);
                RBModel.LF{ii}.FNq{1} = RBModel.LF{ii}.FNq{1}(1:count,:);
                for kk=1:RBModel.LF{ii}.P
                    RBModel.LF{ii}.ANq{kk} = RBModel.LF{ii}.ANq{kk}(1:count,1:count);
                end
                
                RBModel.LF{ii}.N = count;
                disp(['Transformation matrix, ANq matrices and FNq matrix cut at ' num2str(count) ' snapshots.'])
            end
            
            %save ROM
            save([top '/Results/ROM/RBModel.mat'],'RBModel')
        end
        
        function cutMeasurements(varargin)
        %
        %   OrderedModelClass.cutMeasurements(name1,value1,name2,value2,...)
        %
        % Arguments:
        %   num_samples     - the number of samples to cut
        %   num_patterns    - the number of electrode patterns to cut
        %
        %
        
            params = struct();
            for i = 1:2:length(varargin) % work for a list of name-value pairs
                if ischar(varargin{i}) % check if is character
                    params.(varargin{i}) = varargin{i+1}; % override or add parameters to structure.
                end
            end

            data = getenv("ROMEG_DATA");
            if strcmp(data(end),'/'), data = data(1:end-1); setenv("ROMEG_DATA",data); end

            for ii = 1:params.num_samples
                for jj = 1:params.num_patterns
                    OrderedModelClass.changePath(['Result' num2str(ii)])
                    top = getenv("ROMEG_TOP");
                    load([top '/Results/measurements/pattern_' num2str(jj) '.mat'],'Data')
                    Data.p = [];
                    Data.t = [];
                    Data.f = [];
                    save([top '/Results/measurements/pattern_' num2str(jj) '.mat'],'Data')
                    clear Data
                end
                disp(['Done cleaning for sample ' num2str(ii)])
            end
        end
    end
end






















































