function cluster_job_bound()
    
    coef_n=str2double(getenv('SLURM_ARRAY_TASK_ID'));
    top = getenv("ROMEG_TOP");
    fprintf('\n ************************* COEFFICIENT NUMBER  %d \n', coef_n);
    max_snap = str2double(getenv("MAX_SNAP"));
    bound = BoundClass('injection',coef_n,'top',top,'max_snap',max_snap);
    setenv("ROMEG_LOGS",[top '/Results/logs'])
    bound=bound.startLogger(['BOUND_' num2str(coef_n)]);
    bound = bound.makeBound();
    save([top '/Results/bound/bound_' num2str(coef_n) '.mat'],'bound')
    
end