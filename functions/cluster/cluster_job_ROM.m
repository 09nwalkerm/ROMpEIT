function cluster_job_ROM()

electrode=str2double(getenv('SLURM_ARRAY_TASK_ID'));
top = getenv("ROMEG_TOP");
fprintf('\n ************************* ELECTRODE NUMBER  %d \n \n \n', electrode);

ROM=ROMClass('electrode',electrode,'top',top);
setenv("ROMEG_LOGS",[top '/Results/logs'])
ROM=ROM.startLogger(['ROM_' num2str(electrode)]);
ROM=ROM.popFields();
ROM=ROM.buildROM();
ROM.saveLF(electrode)
FOM.logger.info('cluster_job_ROM','Closing MATLAB instance.')
end