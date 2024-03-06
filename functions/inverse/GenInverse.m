function GenInverse(varargin)
%
%   GenInverse(name1,value1,name2,value2...)
%
% Description:
%   Function to run the inverse problem with either the traditional method
%   or using ROM with the RBModel created using the GenRBModel() function
%   (see 'help GenRBModel' for guidance). This function takes the
%   measurements provided with the data_path argument or uses the
%   measurements in the Results/measurements/ folder to estimate the
%   conductivities in the tissues. This function also has the ability to
%   use anisotropic layers if ROM was made with them. 
%
% Arguments:
%   ROM             - Boolean, run ROM inverse problem
%   Cluster         - Boolean, run ROM inverse on the cluster (recommended)
%   TRAD            - Boolean, run the traditional method of estimation
%   verbose         - (beta), Boolean, run in verbose mode
%   simultaneous    - (beta), Boolean, run ROM in simultaneous mode
%   x0              - starting conductivity values for optimization
%   data_path       - path to the measurements
%   model           - path to the model
%   sinks           - 2D array with injection pattern data
%   sinks_path      - path to injection data mat file
%   top             - path to the top of the ROMEG tree
%   current         - the injection current
%   folder          - which Results folder would you like to change into
%   num_samples     - number of samples measurements to run the inverse with
%   sample_num      - the specific sample or samples to run for
%   snaps           - (boolean) would you like the inverse to be run for each
%                     number of snapshots (this will take much longer)
%   fix_conds       - (boolean) fixes the non_active parameters to the middle
%   active_layers   - array of the layers to keep active for inverse
%   sensitivity     - is this a sensitivity analysis?
%   new_sinks       - (boolean) use a new set of sinks called
%                     new_sinks.mat; this means the first set (used for
%                     ROM) must be a 1-1 set where the injection electrode
%                     is the same number as the pattern number and the
%                     extraction electrode is the final electrode.
%   use_sinks       - tell this function that ROM used 'use_sinks'
%   complim         - max number of processes to run in parallel on cluster
%   use_noise       - (boolean) if true, uses measurements with noise
%   noise           - what noise is used in the measurements
%   tag             - (string) how to label the estimates
%

    params = [];
    params_S = [];
    paramslist = [{'ROM'},{'Cluster'},{'TRAD'},{'verbose'},{'simultaneous'}, ...
        {'x0'},{'data_path'},{'model'},{'sinks'},{'sinks_path'},{'top'}, ...
        {'current'},{'folder'},{'num_samples'},{'snaps'},{'fix_conds'}, ...
        {'active_layers'},{'sensitivity'},{'use_sinks'},{'new_sinks'}, ...
        {'complim'},{'use_noise'},{'sample_num'},{'noise'},{'tag'}];

    if ~isempty(varargin)
        for i = 1:2:length(varargin) % work for a list of name-value pairs
            if ischar(varargin{i}) && ~isempty(find(strcmp(paramslist,varargin{i})))
                params = [params varargin(i) varargin(i+1)];
                params_S.(varargin{i}) = varargin{i+1};
            end
        end
    end
    
    if isfield(params_S,'folder')
        OrderedModelClass.changePath(params_S.folder)
    end

    % Get the top of the ROMEG tree
%     if ~isfield(params_S,'top')
%         OrderedModelClass.changePath('Result1')
%         top = getenv("ROMEG_TOP");
%         params = [params {'top'} {top}];
%         
%     end

    if ~isfield(params_S,'model')
        disp("Using example Spherical head model.")
        tree = getenv("ROMEG");
        params=[{'model'} {strcat(tree,'/Models/Spherical/head_model.mat')} params];
    end

    if ~isfield(params_S,'current')
        disp('Using a current of 0.02e-3A') 
        params = [params {'current'} {0.02e-3}];
    end
    
    if ~isfield(params_S,'num_samples') && ~isfield(params_S,'sample_num')
        samples = 1;
    elseif isfield(params_S,'sample_num')
        samples = params_S.sample_num;
    elseif isfield(params_S,'num_samples')
        samples = 1:params_S.num_samples;
    end
    
    if ~isfield(params_S,'active_layers')
        error('Please specify active layer as an array ')
    end

    for i = samples
        if isfield(params_S,'ROM') && params_S.ROM
            
            OrderedModelClass.sensitivityFiles('layers',params_S.active_layers,'sample_num',i,'order','ROM')
            top = getenv("ROMEG_TOP");
            params2 = [params {'top'} {top}];
            
            invROM = InverseROMClass(params2);
            invROM = invROM.checkPaths('type','inverse','num',i);
            invROM.savePrep();
            
            if isfield(params_S,'simultaneous') && params_S.simultaneous
                setenv("num",num2str(1));
                num_injections = 1;
            else
                setenv("num",num2str(invROM.num_patterns));
                num_injections=invROM.num_patterns;
            end
    
            if isfield(params_S,'Cluster') && params_S.Cluster
                
                if isfield(params_S,'complim')
                    setenv("COMPLIM",num2str(params_S.complim))
                    !sbatch --array 1-$num%$COMPLIM -o $ROMEG_TOP/Results/slurm_logs/INVROM_%a_%j.out -e $ROMEG_TOP/Results/slurm_logs/INVROM_%a_%j.err --job-name INVROM $ROMEG/Functions/Cluster/cluster_job.sh inv_ROM
                else
                    !sbatch --array 1-$num -o $ROMEG_TOP/Results/slurm_logs/INVROM_%a_%j.out -e $ROMEG_TOP/Results/slurm_logs/INVROM_%a_%j.err --job-name INVROM $ROMEG/Functions/Cluster/cluster_job.sh inv_ROM
                end
                
                OrderedModelClass.wait('INVROM',20);
            else
                for ii=1:num_injections
                    invROM.runInverse(ii);
                    invROM.saveInv();
                    disp(['Finished pattern ' num2str(ii) ' for sample ' num2str(i)])
                end
            end
            invROM.collect();
        end
    
        if isfield(params_S,'TRAD') && params_S.TRAD
            
            OrderedModelClass.sensitivityFiles('layers',params_S.active_layers,'sample_num',i,'order','TRAD')
    
            invTRAD = InverseTradClass(params);
            invTRAD = invTRAD.checkPaths('type','inverse','num',i);
            invTRAD.savePrep();
            setenv("num",num2str(invTRAD.num_patterns));
    
            if isfield(params_S,'complim')
                setenv("COMPLIM",num2str(params_S.complim))
                !sbatch --array 1-$num%$COMPLIM -o $ROMEG_TOP/Results/slurm_logs/INVTRAD_%a_%j.out -e $ROMEG_TOP/Results/slurm_logs/INVROM_%a_%j.err --job-name INVTRAD $ROMEG/Functions/Cluster/cluster_job.sh inv_TRAD
            else
                !sbatch --array 1-$num -o $ROMEG_TOP/Results/slurm_logs/INVTRAD_%a_%j.out -e $ROMEG_TOP/Results/slurm_logs/INVTRAD_%a_%j.err --job-name INVTRAD $ROMEG/Functions/Cluster/cluster_job.sh inv_TRAD
            end
            
            OrderedModelClass.wait('INVTRAD',30);
            invTRAD.collect();
        end
    end
    
end