%% Ground Test

%% Prelims

%clear

tree = getenv("ROMEG");
%model = [tree '/models/Real/head_model.mat'];
%model2 = '/cubric/data/c1616132/DataDump/ROMEG/Experiments/Real/head_model.mat';
model = '/home/c1616132/Documents/PhD/new_colin/CEM_semi_complete_with_marrow_133_merged.mat';
model2 = '/home/c1616132/Documents/PhD/new_colin/CEM_semi_complete_with_marrow_133.mat';

mu_min = [0.303,0.002,0.013,1.450,0.268,0.092,5];
mu_max = [0.444,0.009,0.043,1.794,0.508,0.177,5];

% For merged skull and spongiform bone
%mu_min = [0.303,0.002,1.450,0.268,0.092,5];
%mu_max = [0.444,0.043,1.794,0.508,0.177,5];

num_samples=4;
num_start=4;

c = [1 2 3 4 5];
%synth = [0.4,0.01,0.01,1.6,0.33,0.12,5];

current_injection = 2e-05;

%% Make measurements
% Type `help GenMeasurements` for examples, help and options

% GenMeasurements('model',model2,'noise',0.82e-6,'mu_min',mu_min,'mu_max',mu_max,...
%     'sample_num',num_start:num_samples,'current',current_injection,...
%     'debug',true,'Cluster',true)

%% Collect Data

%OrderedModelClass.patterns('model',model,'num_sinks',10,'new_sinks',true,'electrodes',1:132)
OrderedModelClass.patterns('model',model,'num_sinks',1,'new_sinks',true,'electrodes',1:132)

data = getenv("ROMEG_DATA");

load([data '/ROM/Results/ROM/sinks.mat'],'sinks')
inv = InverseROMClass('top',[data '/Result' num2str(num_start)],...
    'noise',0.82e-6,'use_noise',true,'debug',true,...
    'ref_sink',133,'simultaneous',true,'new_sinks',true);

inv = inv.collectData();

% for i=1:size(sinks,1)
%    inv.pattern = i;
%    inv = inv.collectData();
% end

load([data '/ROM/Results/ROM/RBModel.mat'],'RBModel')

%% Run the inverse problem

for i=[1 2 5 20]%1:133
    GenInverse('model',model,'ROM',true,'sample_num',num_start:num_samples,...
        'active_layers',c,'use_noise',true,'new_sinks',true,...
        'tag',['1-1_opp_g' num2str(i)],'noise',0.82e-6,'ground',i,'ref_sink',133,...
        'current',current_injection,'u',inv.u,'fix_conds',[5],...
        'RBModel',RBModel,'Display','off')
    disp(['Finished ground electrode ' num2str(i)])
end


%% collect ground test

stdd = zeros(132,5);
for i=1:132
    load([data '/Result' num2str(num_start) '/Results/inverse/ROM/inverse_12345/1-1_opphnr_g' num2str(i) '_estimates.mat'],'estimates')
    indx1 = find(sinks(:,1)==i);
    indx2 = find(sinks(:,2)==i);
    estimates([indx1 indx2],:) = [];
    stdd(i,:) = std(estimates);
end
stdd2 = stdd(:,1)/mean(stdd(:,1)) + stdd(:,2)/mean(stdd(:,2));