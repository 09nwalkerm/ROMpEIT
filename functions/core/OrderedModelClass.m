classdef OrderedModelClass
    properties
        p
        t
        f
        Ind_E           % Tetrahedron number related to electrode
        anis_rad
        anis_tan
        angles = false  % true if model contains theta values (default: false)
        model           % path to model
        theta
        elec_height
        top             % path to top of the ROMEG tree
        sinks           % injection and sink patterns for each electrode
        sinks_path      % path to the sink.mat file
        use_sinks = true % do you want to use the sinks file for ROM/inverse, defualt = true
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
        debug           % (boolean) turn debug mode on for logs with 'true'
        Cluster         % is ROM being run on the cluster
        log_tag
        logger
        close          % see patterns function
        num_weights
        distance
        pre_stiff      % are the stiff mats already in head_model file?
        L
        Aq
    end
    
    %properties (Access = protected)
    %    logger
    %end

    methods
        function obj = OrderedModelClass(varargin)
            obj = obj.processArgs(varargin);
            obj = obj.startLogger();
        end
        
        function obj = getTOP(obj)
            if isempty(obj.top)
                top = getenv("ROMEG_TOP");
                if isempty(top)
                    error('ROMEG_TOP environment variable not set')
                else
                    obj.top = top;
                end
            end
        end
        
        function obj = startLogger(obj)
            
            if ~isempty(obj.log_tag)
                tag = obj.log_tag;
            else
                tag = 'main';
            end
            
            % clear previous logger
            obj.logger = [];
            
            % get env vars and make path
            logs = getenv("ROMEG_LOGS");
            id = getenv("SLURM_JOB_ID");
            logpath=[logs '/' tag '_' id '.log'];
            
            % fetch logger object and load into class
            obj.logger = log4m.getLogger(logpath);
            obj.logger.setCommandWindowLevel(obj.cmdLevel);
            if ~isempty(obj.debug) && obj.debug
                obj.logger.setLogLevel(obj.logger.DEBUG);
                obj.logger.debug('OrderedModelClass.startLogger','Turning on debug mode')
            else
                obj.logger.setLogLevel(obj.logLevel);
            end
            obj.logger.debug('OrderedModelClass.startLogger','Loading logger into Class as obj.logger')
        end
        
        function [M_mu,b_mu] = muAssemble(obj,mu)
            
            if isa(obj,"ROMClass")
                matrix_field = 'ANq';
                source_field = 'FNq';
                obj.logger.trace('muAssemble','Using ANq and FNq matrices in assembly')
            elseif isa(obj,"FOMClass")
                matrix_field = 'Aq';
                source_field = 'Fq';
                obj.logger.debug('muAssemble','Using Aq and Fq matrices in assembly')
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
            
            while length(args)==1 && isa(args,'cell')
                args = args{1};
            end
            
            if ~isempty(args)
                for i = 1:2:length(args) % work for a list of name-value pairs
                    if ischar(args{i}) % check if is character
                        obj.(args{i}) = args{i+1}; % override or add parameters to structure.
                    end
                end
            end
        end

        function obj = processModel(obj)

            obj.logger.info('processModel','Loading Model...')

            
            if ~isempty(obj.pre_stiff) && obj.pre_stiff
                try
                    obj.logger.info('processModel','loading head model with stiffness matrices')
                    load(obj.model,'p','t','L','Aq')
                    obj.p = p; obj.t = t; obj.L = L; obj.Aq = Aq;%obj.Ind_E = Ind_E;
                    obj.logger.info('processModel','p,t,L,Aq from head model loaded')
                catch ME
                    obj.logger.error('processModel',['path given:' obj.model])
                    obj.logger.error('processModel','Cannot load head model, please check path given and that the file contains p,t,L,Aq values')
                    %error('Cannot load head model')
                    rethrow(ME)
                end
                return
            end
            
            if isempty(obj.angles) || (~obj.angles)%isempty(obj.anis_rad) || (~isempty(obj.anis_rad) && (~obj.angles)) %|| ~isfield(obj,'anis_rad')
                try
                    load(obj.model,'p','t','f')
                    obj.p = p; obj.t = t; obj.f = f; %obj.Ind_E = Ind_E;
                    obj.logger.info('processModel','p,t,f from head model loaded')
                catch
                    obj.logger.error('processModel',['path given:' obj.model])
                    obj.logger.error('processModel','Cannot load head model, please check path given and that the file contains p,t,f values')
                    error('Cannot load head model')
                end
            else
                try
                    load(obj.model,'p','t','f','theta')
                    obj.p = p; obj.t = t; obj.f = f; obj.theta = theta;
                    obj.logger.info('processModel','p,t,f,theta from head model loaded with angles loaded')
                catch ME
                    switch ME.identifier
                        case 'MATLAB:load:couldNotReadFile'
                            obj.logger.error('processModel',['path given:' obj.model])
                            obj.logger.error('processModel','Cannot load head model, please check path given')
                            error('Cannot load head model')
                        case 'MATLAB:UndefinedFunction'
                            obj.logger.error('processModel','At least one variable missing, please ensure the head model file contains p,t,f,theta variables.')
                            obj.logger.error('processModel','Or specify the name-value pair angles-false to have it generated.')
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
                            load([obj.top '/Results/ROM/new_sinks.mat'],'sinks')
                            obj.logger.info('loadSinks',['Loading sinks from ' obj.top '/Results/ROM/new_sinks.mat'])
                        catch
                            obj.logger.error('loadSinks','No path or new_sinks.mat file provided and none available in Results folder.')
                            obj.logger.error('loadSinks',"Please make new_sinks.mat using OrderedModelClass.patterns and rerun function.")
                            obj.logger.error('loadSinks',"Alternatively make one manually with custom injection patterns")
                            obj.logger.error('loadSinks','See `help OrderedModelClass.patterns` for guidance.')
                            error('No new_sinks.mat file found')
                        end
                    else
                        try
                            disp(['Loading sinks from ' obj.top '/Results/ROM/sinks.mat'])
                            load([obj.top '/Results/ROM/sinks.mat'],'sinks')
                            obj.logger.info('loadSinks',['Loading sinks from ' obj.top '/Results/ROM/sinks.mat'])
                        catch
                            obj.logger.error('loadSinks','No path or sinks.mat file provided and none available in Results folder.')
                            obj.logger.error('loadSinks',"Please make sinks.mat using OrderedModelClass.patterns and rerun function.")
                            obj.logger.error('loadSinks',"Alternatively make one manually with custom injection patterns")
                            obj.logger.error('loadSinks','See `help OrderedModelClass.patterns` for guidance.')
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
                if strcmp(params.type,'ROM')
                    if ~isfolder([data '/ROM'])
                        disp('Making ROM folder.')
                        mkdir([data '/ROM'])
                        OrderedModelClass.changePath('ROM');
                        %obj.top = [data '/ROM'];
                        OrderedModelClass.setupFiles('ROM',true);
                        obj.logger.debug('checkPaths','Set up ROM folder below $ROMEG_DATA')
                    else
                        %disp('Warning: ROM folder in data path already exists, overwriting.')
                        OrderedModelClass.changePath('ROM')
                        obj.logger.warn('checkPaths','ROM folder already set up, overwriting')
                        %obj.top = [data '/ROM'];
                    end
                elseif strcmp(params.type,'measurement') || strcmp(params.type,'inverse') || strcmp(params.type,'bound')
                    if ~isfolder([data '/Result' num2str(params.num)])
                        disp(['Making Result' num2str(params.num) ' folder.'])
                        mkdir([data '/Result' num2str(params.num)])
                        OrderedModelClass.changePath(['Result' num2str(params.num)])
                        %obj.top = [data '/Result' num2str(params.num)];
                        if isfield(params,'RB_path')
                            OrderedModelClass.setupFiles('ROM',false,'RB_path',params.RB_path)
                        else
                            OrderedModelClass.setupFiles('ROM',false)
                        end
                        obj.logger.debug('checkPaths',['Set up file system in Result' num2str(params.num) 'below $ROMEG_DATA'])
                    else
                        obj.logger.warn('checkPaths',['Result' num2str(params.num) ' folder in data path already exists, overwriting.'])
                        OrderedModelClass.changePath(['Result' num2str(params.num)])
                        %obj.top = [data '/Result' num2str(params.num)];

                    end
                elseif strcmp(params.type,'eeg') || strcmp(params.type,'LF')
                    if ~isfolder([data '/Result' num2str(params.num)])
                        disp(['Making Result' num2str(params.num) ' folder.'])
                        mkdir([data '/Result' num2str(params.num)])
                        OrderedModelClass.changePath(['Result' num2str(params.num)])
                        if strcmp(params.type,'LF'); params.num_dipoles=0; end
                        OrderedModelClass.EEGFiles('sample_num',params.num,'num_dipoles',params.num_dipoles);
                    else
                        obj.logger.warn('checkPaths',['Result' num2str(params.num) ' folder in data path already exists, overwriting.'])
                        OrderedModelClass.changePath(['Result' num2str(params.num)])
                        if strcmp(params.type,'LF'); params.num_dipoles=0; end
                        OrderedModelClass.EEGFiles('sample_num',params.num,'num_dipoles',params.num_dipoles);
                        %obj.top = [data '/Result' num2str(params.num)];
                    end
                else
                    OrderedModelClass.changePath('ROM');
                end
                obj.top = getenv("ROMEG_TOP");
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
        %       model: (essential) path to the head model
        %       num_sinks: how many current sinks for each
        %                     injection would you like.
        %       elec_height: (optional) a given minimum height the sink
        %                     electrodes should be on the head model.
        %       electrodes: list of electrodes for injection (defaults to
        %                     all)
        %       top: (optional if ROMEG_TOP env var has been set)
        %                     path to the top of the ROMEG tree
        %       out: specific list of extraction electrodes for
        %                     injection patterns.
        %       new_sinks: are these new sinks being made for inverse
        %       close: select the closest electrodes. Default is furthest
        %                     away
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

            obj = OrderedModelClass(varargin);
            obj = obj.checkPaths('type','ROM');
            obj = obj.processModel();
            f = obj.f; p = obj.p;

            obj.logger.info('patterns','Generating sink patterns...')
            
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

                    for jj = 1:size(elec_sinks,1)
                        lengths(jj,:) = elec_sinks(jj,:) - pos(:,:);
                        if ~isempty(obj.elec_height)
                            if (elec_sinks(jj,3) < obj.elec_height), lengths(jj,:) = []; end
                        end
                    end

                    norms = vecnorm(lengths');
                    
                    if ~isempty(obj.close) && obj.close
                        [~,ind] = mink(norms,obj.num_sinks+1);
                        ind = ind(2:end);
                    else
                        [~,ind] = maxk(norms,obj.num_sinks);
                    end

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
                
            obj.logger.info('patterns','Saved sink patterns to /Results/ROM folder')
        end
        
        function weights = makeWeights(varargin)
        %
        %   OrderedModelClass.makeWeights(name1,value1,name2,value2...)
        %
        %   Description:
        %       A script to make the weights for the measurements on each
        %       eelctrode for each injection pattern in the inverse
        %       problem.
        % 
        %   Arguments:
        %       model: (essential) path to the head model
        %       num_weights: how many current sinks for each
        %                     injection would you like.
        %       distance: should the weights be a function of distance from
        %                 the injection electrode.
        %       sinks: matrix containing all the injection patterns
            
            obj = OrderedModelClass(varargin);
            %obj = obj.checkPaths('type','ROM');
            obj = obj.processModel();
            f = obj.f; p = obj.p;
            
            %load([obj.top '/Results/ROM/new_sinks.mat'],'sinks')
            
            elec_centers = zeros(length(unique(f(:,4)))-1,3);
            
            for ii = 1:length(elec_centers)
                nodes = f(f(:,4)==ii,1:3);
                nodes = unique(nodes);
                elec_centers(ii,:) = [mean(p(nodes,1)) mean(p(nodes,2)) mean(p(nodes,3))];
            end

            weights = zeros(size(obj.sinks,1),length(unique(f(:,4)))-1-size(obj.sinks,2));
            
            for ii = 1:size(obj.sinks,1)
                weight = zeros(1,length(unique(f(:,4)))-1);
                pos = elec_centers(obj.sinks(ii,1),:);

                elec_sinks = elec_centers;
                lengths = zeros(size(elec_sinks,1),3);

                for jj = 1:size(elec_sinks,1)
                    lengths(jj,:) = elec_sinks(jj,:) - pos(:,:);
                    if ~isempty(obj.elec_height)
                        if (elec_sinks(jj,3) < obj.elec_height), lengths(jj,:) = []; end
                    end
                end

                norms = vecnorm(lengths');

                if ~isempty(obj.num_weights)
                    norms(obj.sinks(ii,1)) = NaN;
                    norms(obj.sinks(ii,2:end)) = NaN;
                    [~,ind] = mink(norms,obj.num_weights);
                    weight(ind) = 1;
                end
                
                if ~isempty(obj.distance) && obj.distance
                    max_distance = max(norms);
                    weights(ii,:) = max_distance-norms;
                end
                
                weight = normalize(weight,"norm",1);
                weight(obj.sinks(ii,:)) = [];
                weights(ii,:) = weight;
            end
        end
        
        function wait(NAME,PAUSE,varargin)
            pause(5)
            if nargin > 2
                JOB=varargin{1};
            else
                [~,cmdout] = system(['squeue -u $USER | grep -m 1 ' NAME ' | awk ''{split($1,a,"_"); print a[1]; exit}'' ']);
                cmdout2 = cmdout(1:end-1);
                JOB = cmdout2;
            end
            fprintf('\nWaiting for slurm job to finish')
            while 1
                fprintf('.')
                [~,cmdout] = system(['squeue --job ' JOB ' | grep -m 1 ' NAME ]);
                if isempty(cmdout)
                    [~,cmdout2] = system(['scontrol show job ' JOB ' |grep -B 3 ''JobState=FAILED'' ']);
                    if isempty(cmdout2)
                        break
                    end
                elseif strcmp(cmdout(1:end-1),'slurm_load_jobs error: Invalid job id specified')
                    disp(cmdout)
                    break
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
            data = getenv("ROMEG_DATA");
            if isempty(data)
                error('Must set ROMEG_DATA. Please source set_env.sh.')
            end
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
                !ln -s $ROMEG_DATA/ROM/Results/ROM $ROMEG_TOP/Results/ROM
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
        
    end
end






















































