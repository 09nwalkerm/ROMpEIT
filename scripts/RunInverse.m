%% RunInverse

% A script to run the inverse problem given an RBModel exists in the path
% specified by the evironment variable $ROMEG_DATA

% Optional: create synthetic measurements before running the inverse
% problem.

%% Prelims

tree = getenv("ROMEG");
model = [tree '/Models/Real/head_model.mat'];

data = getenv("ROMEG_DATA");
jobid = getenv("SLURM_JOB_ID");
logger = log4m.getLogger([data '/logs/' jobid '.log']);
logger.setLogLevel(logger.INFO); % set to logger.OFF for only slurm log output
logger.setCommandWindowLevel(logger.INFO); % set to logger.OFF for only log file input

mu_min = [0.303,0.002,0.013,1.450,0.268,0.092,5];
mu_max = [0.444,0.009,0.043,1.794,0.508,0.177,5];

% Anisotropic values in the skull
%mu_min = [0.303,0.002,0.002,1.450,0.268,0.092,5];
%mu_max = [0.444,0.043,0.043,1.794,0.508,0.177,5];

% For merged skull and spongiform bone
%mu_min = [0.303,0.002,1.450,0.268,0.092,5];
%mu_max = [0.444,0.043,1.794,0.508,0.177,5];

num_samples=70;
num_start=60;

c = [1 2 3 4 5 6];
%synth = [0.4,0.01,0.01,1.6,0.33,0.12,5];

%% Make measurements
% Type `help GenMeasurements` for examples, help and options

% GenMeasurements('model',model,'noise',0.82e-6,'mu_min',mu_min,'mu_max',mu_max,...
%     'sample_num',num_start:num_samples)

%% Run the inverse problem

model = [tree '/Models/Real/head_model_theta.mat'];

%GenInverse('model',model,'ROM',true,'Cluster',true,'sample_num',num_start:num_samples,...
%    'active_layers',c,'use_sinks',true,'use_noise',true,'fix_conds',true)

OrderedModelClass.patterns('model',model,'out',82,'new_sinks',true,'electrodes',1:132)

GenInverse('model',model,'ROM',true,'Cluster',true,'sample_num',num_start:num_samples,...
    'active_layers',c,'new_sinks',true,'use_noise',true,'fix_conds',true,...
    'tag','1-1-fixed82','noise',0.82e-6)

% GenInverse('model',model,'ROM',true,'Cluster',true,'sample_num',num_start:num_samples,...
%     'active_layers',c,'use_sinks',true,'use_noise',true,'fix_conds',true,'simultaneous',true,...
%     'tag','1-1_sim','noise',0.82e-6)
% 
% GenInverse('model',model,'ROM',true,'Cluster',true,'sample_num',num_start:num_samples,...
%     'active_layers',c,'new_sinks',true,'use_noise',true,'fix_conds',true,...
%     'tag','1-20','noise',0.82e-6)
% 
% GenInverse('model',model,'ROM',true,'Cluster',true,'sample_num',num_start:num_samples,...
%     'active_layers',c,'new_sinks',true,'use_noise',true,'fix_conds',true,'simultaneous',true,...
%     'tag','1-20_sim','noise',0.82e-6)
% 
% OrderedModelClass.patterns('model',model,'num_sinks',20,'new_sinks',true,'electrodes',1:132,'elec_height',0.05)
% 
% GenInverse('model',model,'ROM',true,'Cluster',true,'sample_num',num_start:num_samples,...
%     'active_layers',c,'new_sinks',true,'use_noise',true,'fix_conds',true,...
%     'tag','1-20_height','noise',0.82e-6)