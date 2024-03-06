%% Script for generating cond map and comparing EEG forward problem

% This script runs the inverse problem to find the conductivities for
% a head model given a set of measurements using the anisotropically
% trained ROM. A map is then found by interpolating the estimated
% tangential and radial conductivities. This model is then used for the EEG
% forward problem given a source of activity in the brain and compared with
% the forward problem using a model with spongiform bone. The error between
% the two is then compared with the error between the spongiform model and
% a single skull layer isotropic model.

tree = getenv("ROMEG");
data = getenv("ROMEG_DATA");
%load([tree '/Models/Real/head_model.mat'])
jobid = getenv("SLURM_JOB_ID");

logger = log4m.getLogger([data '/logs/' jobid '.log']);
logger.setLogLevel(logger.TRACE); % set to logger.OFF for only slurm log output
logger.setCommandWindowLevel(logger.INFO); % set to logger.OFF for only log file input

%% Prelims
sample_num = 14:20;

% Source position
%pos = [0.06,0.14,0.135;0.06,0.08,0.135;0.1,0.08,0.135];
%q=[1,0,0; 1,0,0; 1,0,0]*1e-8; % Tangentially-oriented dipolar source [A.m]
[pos,q] = sources(100);


for i=sample_num
    model = [tree '/Models/Real/head_model.mat'];
    load(['/cubric/scratch/c1616132/ROMEG_R12345_CUSTOM1/Result' num2str(i) '/Results/measurements/prep.mat'],'Data')
    conds = Data.synth_cond(1:6);
    tag = 'FULL';
    sig_inf1 = Data.synth_cond(5);
    
    GenEEG('model',model,'pos',pos,'q',q,'tag',tag,'sample_num',...
        i,'conds',conds,...
        'sig_inf',sig_inf1)

    model2 = [tree '/Models/Real/head_model_theta.mat'];
    load(['/cubric/scratch/c1616132/ROMEG_R12345_CUSTOM1/Result' num2str(i) '/Results/inverse/ROM/inverse_12345/estimate.mat'],'estimate')
    tag = 'ISH';
    sig_inf = estimate(4);
    
    GenEEG('model',model2,'pos',pos,'q',q,'tag',tag,'sample_num',...
        i,'conds',estimate,...
        'sig_inf',sig_inf)

    load(['/cubric/scratch/c1616132/ROMEG_R12345_CUSTOM1/Result' num2str(i) '/Results/inverse/ROM/inverse_12345/estimates.mat'],'estimates')
    tag = 'ISI';
    sig_inf = mean(estimates(:,4));
    
    GenEEG('model',model2,'pos',pos,'q',q,'tag',tag,'sample_num',...
        i,'conds',estimates,...
        'sig_inf',sig_inf)
    
    load(['/cubric/scratch/c1616132/ROMEG_R12345_CUSTOM2/Result' num2str(i) '/Results/inverse/ROM/inverse_12345/estimates.mat'],'estimates')
    tag = 'IMI';
    sig_inf = mean(estimates(:,4));
    
    GenEEG('model',model2,'pos',pos,'q',q,'tag',tag,'sample_num',...
        i,'conds',estimates,...
        'sig_inf',sig_inf)
    
    load(['/cubric/scratch/c1616132/ROMEG_R12A3A456_CUSTOM5/Result' num2str(i) '/Results/inverse/ROM/inverse_123456/estimates.mat'],'estimates')
    tag = 'AMI';
    sig_inf = mean(estimates(:,5));
    
    GenEEG('model',model2,'pos',pos,'q',q,'tag',tag,'sample_num',...
        i,'conds',estimates,'angles',true,...
        'sig_inf',sig_inf)

end
%% Load FP solutions and plot error

% data = getenv("ROMEG_DATA");
% sample_num = 1;
% dp_num = 2;
% 
% load([data '/Result' num2str(sample_num) '/Results/EEG_FP/DP_' num2str(dp_num) '/full.mat'],'u_n')
% % u_n_full = u_n(1:end-133,:)
% u_n_full = u_n;
% clear u_n
% 
% load([data '/Result' num2str(sample_num) '/Results/EEG_FP/DP_' num2str(dp_num) '/aniso.mat'],'u_n')
% % u_n_aniso = u_n(1:end-133,:);
% u_n_aniso = u_n;
% clear u_n
% 
% load([data '/Result' num2str(sample_num) '/Results/EEG_FP/DP_' num2str(dp_num) '/iso.mat'],'u_n')
% % u_n_iso = u_n(1:end-133,:);
% u_n_iso = u_n;
% clear u_n
% 
% error_iso = abs((u_n_full - u_n_iso)./u_n_full);
% error_aniso = abs((u_n_full - u_n_aniso)./u_n_full);
% 
% norm(error_iso)
% norm(error_aniso)
% 
% %error_scalp_aniso = error_aniso(f(:,1));
% %error_scalp_iso = error_iso(f(:,1));
% 
% %mean(error_scalp_iso)
% %mean(error_scalp_aniso)
% 
% plotting = PlottingClass();
% 
% plotting.plotHeadModel('model',model2,'tissue',4,'fill',error_iso,'cmap','hot')
% plotting.plotHeadModel('model',model2,'tissue',4,'fill',error_aniso,'cmap','hot')




