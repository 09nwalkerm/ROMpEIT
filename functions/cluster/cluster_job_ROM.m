function cluster_job_ROM()

electrode=str2double(getenv('SLURM_ARRAY_TASK_ID'));
top = getenv("ROMEG_TOP");
fprintf('\n ************************* ELECTRODE NUMBER  %d \n \n \n', electrode);

setenv("ROMEG_LOGS",[top '/Results/logs'])
ROM=ROMClass('electrode',electrode,'top',top,'log_tag',['ROM_' num2str(electrode)]);
ROM=ROM.buildROM();
ROM.saveLF(electrode)
ROM.logger.info('cluster_job_ROM','Closing MATLAB instance.')
end