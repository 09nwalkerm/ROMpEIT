%% RunInverse

% A script to run the inverse problem given an RBModel exists in the path
% specified by the evironment variable $ROMEG_DATA

% Optional: create synthetic measurements before running the inverse
% problem.

%% Prelims

tree = getenv("ROMEG");
model = [tree '/models/Real/head_model.mat'];

%mu_min = [0.303,0.002,0.013,1.450,0.268,0.092,5];
%mu_max = [0.444,0.009,0.043,1.794,0.508,0.177,5];

% Anisotropic values in the skull
%mu_min = [0.303,0.002,0.002,1.450,0.268,0.092,5];
%mu_max = [0.444,0.043,0.043,1.794,0.508,0.177,5];

% For merged skull and spongiform bone
%mu_min = [0.303,0.002,1.450,0.268,0.092,5];
%mu_max = [0.444,0.043,1.794,0.508,0.177,5];

num_samples=100;
num_start=1;

c = [1 2 3 4 5 6];
%synth = [0.4,0.01,0.01,1.6,0.33,0.12,5];

%% Make measurements
% Type `help GenMeasurements` for examples, help and options

% GenMeasurements('model',model,'noise',0.82e-6,'mu_min',mu_min,'mu_max',mu_max,...
%     'sample_num',num_start:num_samples)

%% Run the inverse problem

%OrderedModelClass.patterns('model',model,'num_sinks',1,'new_sinks',true,'electrodes',5:132,'close',true)

% GenInverse('model',model,'ROM',true,'Cluster',true,'sample_num',num_start:num_samples,...
%     'active_layers',c,'use_noise',true,...
%     'tag','1-133-known','noise',0.82e-6,'debug',true)

GenInverse('model',model,'ROM',true,'Cluster',true,'sample_num',num_start:num_samples,...
    'active_layers',c,'use_noise',true,'use_sinks',true,'fix_conds',true,...
    'tag','','noise',0.82e-6,'debug',true,'ref_sink',133,'snaps',true)
