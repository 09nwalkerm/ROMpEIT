function GenBound(varargin)
%
%   GenBound(name1,value1,name2,value2...)
%
% Description:
%   Function to generate the necessary data for the bound plot and save it
%   in the 'bound' folder in each Result directory. This function takes the
%   max and mean error bounds calculated during the ROM process and saves
%   those in the BoundClass file along with the reduced basis solution of
%   the forward problem for the same conductivity as the measurements in
%   the same folder.
%
% Arguments:
%   num_samples: number of conductivity samples to use
%   sample_num: array of sample numbers
%   max_snap: maximum number of snapshots in RBModel to use for
%                     calculations
%   complim: max number of processes to run in parallel on cluster
%   debug: (boolean) turn debug mode on


    params = [];
    params_S = struct();
    paramslist = [{'num_samples'},{'max_snap'},{'complim'},{'sample_num'},...
        {'debug'}];

    if ~isempty(varargin)
        for i = 1:2:length(varargin) % work for a list of name-value pairs
            if ischar(varargin{i}) && ~isempty(find(strcmp(paramslist,varargin{i})))
                params = [params varargin(i) varargin(i+1)];
                params_S.(varargin{i}) = varargin{i+1};
            end
        end
    end
    
    if ~isfield(params_S,'num_samples') && ~isfield(params_S,'sample_num')
        samples = 1;
    elseif isfield(params_S,'sample_num')
        samples = params_S.sample_num;
    elseif isfield(params_S,'num_samples')
        samples = 1:params_S.num_samples;
    end

    for i = samples
        
        OMC = OrderedModelClass();
        OMC = OMC.checkPaths('type','bound','num',i);
        OMC = OMC.startLogger('BOUND');
        OMC = OMC.loadSinks();
        setenv("num",num2str(OMC.num_patterns));
        setenv("MAX_SNAP",num2str(params_S.max_snap))
        
        if isfield(params_S,'complim')
            setenv("COMPLIM",num2str(params_S.complim))
            !sbatch --array 1-$num%$COMPLIM -o $ROMEG_TOP/Results/slurm_logs/BOUND_%a_%j.out --job-name BOUND $ROMEG/Functions/Cluster/cluster_job.sh bound
        else
            !sbatch --array 1-$num -o $ROMEG_TOP/Results/slurm_logs/BOUND_%a_%j.out --job-name BOUND $ROMEG/Functions/Cluster/cluster_job.sh bound
        end
        
        OrderedModelClass.wait('BOUND',20);
        
    end

    OrderedModelClass.changePath('ROM'); top = getenv("ROMEG_TOP");
    BC = BoundClass('top',top,'num_samples',samples(end),'max_snap',params_S.max_snap);
    BC = BC.collectBound();
    BC = BC.processBound();
    BC = BC.saveBound();
end