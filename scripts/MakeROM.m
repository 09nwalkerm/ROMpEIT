%% Generate a RBModel

% Script to determine which layers are best to estimate by the reduction in
% overall estimation error.

%% Prelims
tree = getenv("ROMEG");
%model = [tree '/models/Spherical/head_model.mat']; % path to model to use
model = '/home/c1616132/Documents/PhD/colin_sources/head_model_CEM.mat';

mu_min = [0.1359,0.0008,0.0010,1.3881,0.0598,0.0652,5];
mu_max = [0.6196,0.0131,0.0430,1.7945,0.7391,0.2283,5];

current_injection = 2.0e-05;

%% Set up files

OMC = OrderedModelClass();
OMC = OMC.checkPaths('type','ROM');

%% Create sink file
% Type `help OrderedModelClass.patterns` for help, options and examples

OrderedModelClass.patterns('model',model,'electrodes',1:165,'out',165)

%% Run functions to make RB model
% Type `help GenRBModel` for help, options and examples

GenRBModel('model',model,'mu_min',mu_min,'mu_max',mu_max,'tolGREEDY',5e-10,'Nmax',100,...
    'use_sinks',true,'nic',150,'debug',true,...
    'Cluster',true,'current',current_injection)