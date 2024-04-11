function cluster_job_bound()
    
    coef_n=str2double(getenv('SLURM_ARRAY_TASK_ID'));
    top = getenv("ROMEG_TOP");
    fprintf('\n ************************* COEFFICIENT NUMBER  %d \n', coef_n);
    max_snap = str2double(getenv("MAX_SNAP"));
    debug = getenv('DEBUG');
    if strcmp(debug,'true')
        debug = true;
    else
        debug = false;
    end
    setenv("ROMEG_LOGS",[top '/Results/logs'])
    bound = BoundClass('pattern',coef_n,'top',top,'max_snap',max_snap,'log_tag',['BOUND_' num2str(coef_n)],'debug',debug);
    bound.logger.info('cluster_job_bound','Set up BoundClass, now running')
    bound = bound.makeBound();
    save([top '/Results/bound/bound_' num2str(coef_n) '.mat'],'bound')
    %bound.logger.info('cluster_job_bound',['Saved result to ' top '/Results/bound/bound_' num2str(coef_n) '.mat'])
end