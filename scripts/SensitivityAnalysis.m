%% Sensitivity Analysis

% Script to determine which layers are best to estimate by the reduction in
% overall estimation error.

%% Prelims
tree = getenv("ROMEG");
%model = [tree '/Models/Real/head_model.mat']; % path to model to use
model = [tree '/Models/Real/head_model_theta.mat'];

% Conds for spherical head model, comment out as necessary
%mu_max=[0.66,0.06,2.3,1.00,1.00,5]; % maximum conductivities
%mu_min=[0.15,0.001,1.1,0.05,0.05,5]; % minimum conductivities

% Conds for real head model, comment out as necessary
%mu_min=[0.15,0.001,0.001,1.1,0.05,0.05,5];
%mu_max=[0.66,0.006,0.060,2.3,1.00,0.65,5];

mu_min = [0.303,0.002,0.002,1.450,0.268,0.092,5];
mu_max = [0.444,0.043,0.043,1.794,0.508,0.177,5];

num_samples = 50;
num_start = 5;

layer_names = {'Scalp','Skull Tangential','Skull Radial','CSF','Grey Matter','White Matter'};

c = [1 2 3 4 5 6];

%% Run function to make synthetic measurements
% Type `help GenMeasurements` for help, options and examples

GenMeasurements('model',model,'sample_num',num_start:num_samples,'mu_min',mu_min,'mu_max',mu_max,...
    'anis_tan',2,'anis_rad',3,'angles',true,'ratio',5,'noise',0.82e-6)

%% Running the inverse problem
% Type `help GenInverse` for help, options and examples

GenInverse('model',model,'ROM',true,'Cluster',true,'sample_num',num_start:num_samples,...
    'active_layers',c,'sensitivity',true,'use_sinks',true,'use_noise',true)

sense = SenseClass('layers',c,'layer_names',layer_names,'plot','Box');
sense = sense.readResults(c);
sense = sense.processResults(c);
sense = sense.saveSense();

% Type `GenPlot('sense',true)` for plot
