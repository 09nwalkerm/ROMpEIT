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
%   ROM: Boolean, run ROM inverse problem
%   Cluster: Boolean, run ROM inverse on the cluster (recommended)
%   TRAD: Boolean, run the traditional method of estimation
%   simultaneous: (beta), Boolean, run ROM in simultaneous mode
%   x0: starting conductivity values for optimization
%   data_path: path to the measurements
%   model: path to the model
%   sinks: 2D array with injection pattern data
%   sinks_path: path to injection data mat file
%   top: path to the top of the ROMEG tree
%   current: the injection current
%   num_samples: number of samples measurements to run the inverse with
%   sample_num: the specific sample or samples to run for
%   snaps: (boolean) would you like the inverse to be run for each
%                     number of snapshots (this will take much longer)
%   fix_conds: (boolean) fixes the non_active parameters to the middle
%   active_layers: array of the layers to keep active for inverse
%   sensitivity: is this a sensitivity analysis?
%   new_sinks: (boolean) use a new set of sinks called
%                     new_sinks.mat; this means the first set (used for
%                     ROM) must be a 1-1 set where the injection electrode
%                     is the same number as the pattern number and the
%                     extraction electrode is the final electrode.
%   use_sinks: tell this function that ROM used 'use_sinks'
%   complim: max number of processes to run in parallel on cluster
%   use_noise: (boolean) if true, uses measurements with noise
%   noise: what noise is used in the measurements
%   tag: (string) how to label the estimates
%   debug: (boolean) turn debug mode on
%   weighted: use weights on the measurements
%   ref_sink: used when ROM is trained with a common sink (advised), and
%       then new sinks are used in the inverse problem
%   omit_layers: when using synthetic measurements with different number of
%       layers. Specify the layer(s) not trained in ROM model.
%

    params = [];
    params_S = [];
    paramslist = [{'ROM'},{'Cluster'},{'TRAD'},{'simultaneous'}, ...
        {'x0'},{'data_path'},{'model'},{'sinks'},{'sinks_path'},{'top'}, ...
        {'current'},{'num_samples'},{'snaps'},{'fix_conds'}, ...
        {'active_layers'},{'sensitivity'},{'use_sinks'},{'new_sinks'}, ...
        {'complim'},{'use_noise'},{'sample_num'},{'noise'},{'tag'},...
        {'debug'},{'ref_sink'},{'weighted'},{'omit_layers'},{'iter'}];

    if ~isempty(varargin)
        for i = 1:2:length(varargin) % work for a list of name-value pairs
            if ischar(varargin{i}) && ~isempty(find(strcmp(paramslist,varargin{i})))
                params = [params varargin(i) varargin(i+1)];
                params_S.(varargin{i}) = varargin{i+1};
            end
        end
    end

    if ~isfield(params_S,'model')
        error('Please specify a path to a head model')
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
        error('Please specify active layers as an array ')
    end
    
    [invROM,invTRAD] = prep(params_S,params,samples);
    [invROM,invTRAD] = run(params_S,samples,invROM,invTRAD);
    
    if isfield(params_S,'iter')
        for i=samples
            k = dbscan(invROM{i}.estimates(:,1),0.0005,2);
            invROM{i}.sinks(k~=1,:) = [];
            invROM{i}.logger.info('GenInverse',['Cut sample ' num2str(i) ' patterns file down to ' num2str(size(invROM{i}.sinks,1)) ' patterns'])
            invROM{i}.num_patterns = size(invROM{i}.sinks,1);
            invROM{i}.weights(k~=1,:) = [];
            invROM{i}.simultaneous = true;
            invROM{i}.savePrep();
        end
        params_S.simultaneous = true;
        [~,~] = run(params_S,samples,invROM,invTRAD);
    end
end

function [invROM,invTRAD] = prep(params_S,params,samples)
    invROM = [];
    invTRAD = [];
    for i = samples
        if isfield(params_S,'ROM') && params_S.ROM
            
            OrderedModelClass.sensitivityFiles('layers',params_S.active_layers,'sample_num',i,'order','ROM')
            
            inv = InverseROMClass(params);
            inv = inv.checkPaths('type','inverse','num',i);
            inv.savePrep();
            invROM{i} = inv;
        end
        
        if isfield(params_S,'TRAD') && params_S.TRAD
            
            OrderedModelClass.sensitivityFiles('layers',params_S.active_layers,'sample_num',i,'order','TRAD')
    
            inv = InverseTradClass(params);
            inv = inv.checkPaths('type','inverse','num',i);
            inv.savePrep();
            invTRAD{i} = inv;
        end
    end
end

function [invROM,invTRAD] = run(params_S,samples,invROM,invTRAD)
    for i = samples
        setenv("ROMEG_TOP",invROM{i}.top)
        if isfield(params_S,'ROM') && params_S.ROM
            
            if isfield(params_S,'simultaneous') && params_S.simultaneous
                setenv("num",num2str(1));
                num_injections = 1;
            else
                setenv("num",num2str(invROM{i}.num_patterns));
                num_injections = invROM{i}.num_patterns;
                invROM{i}.logger.info('run',['setting num_injections to ' num2str(num_injections)])
            end
            
            if isfield(params_S,'Cluster') && params_S.Cluster
                
                if isfield(params_S,'complim')
                    setenv("COMPLIM",num2str(params_S.complim))
                    !sbatch --array 1-$num%$COMPLIM -o $ROMEG_TOP/Results/slurm_logs/INVROM_%a_%j.out -e $ROMEG_TOP/Results/slurm_logs/INVROM_%a_%j.err --job-name INVROM $ROMEG/functions/cluster/cluster_job.sh inv_ROM
                else
                    !sbatch --array 1-$num -o $ROMEG_TOP/Results/slurm_logs/INVROM_%a_%j.out -e $ROMEG_TOP/Results/slurm_logs/INVROM_%a_%j.err --job-name INVROM $ROMEG/functions/cluster/cluster_job.sh inv_ROM
                end
                
                if ~isfield(params_S,'simultaneous')
                    OrderedModelClass.wait('INVROM',20);
                else
                    pause(3)
                end
            else
                for ii=1:num_injections
                    invROM{i}.runInverse(ii);
                    invROM{i}.saveInv();
                    invROM{i}.logger.info('GenInverse',['Finished ROM pattern ' num2str(ii) ' for sample ' num2str(i)])
                end
            end
            try
                if ~isfield(params_S,'simultaneous')
                    invROM{i} = invROM{i}.collect();
                end
            catch ME
                invROM{i}.logger.error('GenInverse','Collection of estimates failed')
                rethrow(ME)
            end
        end
    
        if isfield(params_S,'TRAD') && params_S.TRAD
            if isfield(params_S,'Cluster') && params_S.Cluster
                
                setenv("num",num2str(invTRAD.num_patterns));
                if isfield(params_S,'complim')
                    setenv("COMPLIM",num2str(params_S.complim))
                    !sbatch --array 1-$num%$COMPLIM -o $ROMEG_TOP/Results/slurm_logs/INVTRAD_%a_%j.out -e $ROMEG_TOP/Results/slurm_logs/INVROM_%a_%j.err --job-name INVTRAD $ROMEG/functions/cluster/cluster_job.sh inv_TRAD
                else
                    !sbatch --array 1-$num -o $ROMEG_TOP/Results/slurm_logs/INVTRAD_%a_%j.out -e $ROMEG_TOP/Results/slurm_logs/INVTRAD_%a_%j.err --job-name INVTRAD $ROMEG/functions/cluster/cluster_job.sh inv_TRAD
                end

                OrderedModelClass.wait('INVTRAD',30);
            else
                for ii=1:invTRAD{i}.num_patterns
                    invTRAD{i}.runInverse(ii);
                    invTRAD{i}.saveInv();
                    disp(['Finished TRAD pattern ' num2str(ii) ' for sample ' num2str(i)])
                end
            end
            try
                invTRAD{i} = invTRAD{i}.collect();
            catch ME
                invTRAD{i}.logger.error('GenInverse','Collection of estimates failed')
                rethrow(ME)
            end
        end
    end
    
    if isfield(params_S,'simultaneous') && params_S.simultaneous
        OrderedModelClass.wait('INVROM',20);
        for s=samples
            invROM{s} = invROM{s}.collect();
        end
    end
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    