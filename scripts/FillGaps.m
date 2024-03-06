%% Filling the gaps caused by faulty nodes
%
% This script will find and re-run any gaps in the INVTRAD results that
% have come from faulty nodes. The --no-requeue option is in the sbatch
% scripts to stop the work flow being disrupted.

tree = getenv("ROMEG");
data = getenv("ROMEG_DATA");
%load([tree '/Models/Real/head_model.mat'])
jobid = getenv("SLURM_JOB_ID");

logger = log4m.getLogger([data '/logs/' jobid '.log']);
logger.setLogLevel(logger.INFO); % set to logger.OFF for only slurm log output
logger.setCommandWindowLevel(logger.INFO); % set to logger.OFF for only log file input

num_end=15;
num_start=11;
num_elecs = 132;
inv_folder = 'inverse_123';

samples=num_start:num_end;

for ii=1:length(samples)
    folder = ['Result' num2str(samples(ii))];
    OrderedModelClass.changePath(folder)
    top = getenv("ROMEG_TOP");
    load([top '/Results/inverse/TRAD/' inv_folder '/estimate.mat'],'estimates','histories')
    for jj = 1:num_elecs
        if (estimates(jj,:) == [0 0 0])%~isfile([top '/Results/inverse/TRAD/' inv_folder '/inv_' num2str(jj) '.mat'])
            setenv("num",num2str(jj));
            !sbatch --array $num -o $ROMEG_TOP/Results/slurm_logs/INVTRAD_%a_%j.out -e $ROMEG_TOP/Results/slurm_logs/INVTRAD_%a_%j.err --job-name INVTRAD $ROMEG/Functions/Cluster/cluster_job.sh inv_TRAD
        end
    end
    OrderedModelClass.wait('INVTRAD',30);
    pause(10)
    %load([top '/Results/inverse/TRAD/' inv_folder '/estimate.mat'],'estimates','histories')
    for jj=1:num_elecs
        if isfile([top '/Results/inverse/TRAD/' inv_folder '/inv_' num2str(jj) '.mat'])
            load([top '/Results/inverse/TRAD/' inv_folder '/inv_' num2str(jj) '.mat'],'inv')
            estimates(jj) = inv.estimate;
            histories{jj} = inv.history;
            %delete([top '/Results/inverse/TRAD/' inv_folder '/inv_' num2str(jj) '.mat'])
        end
    end
    save([top '/Results/inverse/TRAD/' inv_folder '/estimate2.mat'],'estimates','histories')
end


