%% Data check

data = getenv("ROMEG_DATA");

folder = 106;
ref = 133;

load([data '/Result' num2str(folder) '/Results/inverse/ROM/inverse_123456/inc_sutures_estimates.mat'])
load([data '/ROM/Results/ROM/RBModel.mat'])
load([data '/ROM/Results/ROM/new_sinks.mat'])

for rr = 53%1:size(sinks,1)
    
    % Conductivities
    cond=[estimates(rr,1);estimates(rr,2);estimates(rr,3);estimates(rr,4);estimates(rr,5);estimates(rr,6)]; 

    % Contact impedances
    L=133;

    % Current injection pattern
    ubi1=sinks(rr,1); ubi2=sinks(rr,2);
    
    % solution for injection
    n_mu=7; % number of parameters
    M_mu=RBModel.LF{ubi1}.ANq{n_mu}/5; % Z

    for kk=1:n_mu-1
        M_mu=M_mu+cond(kk)*RBModel.LF{ubi1}.ANq{kk}; % Stiffness
    end

    zN1=M_mu\RBModel.LF{ubi1}.FNq{1};
    zNh1 = RBModel.LF{ubi1}.V*zN1;

    clear M_mu

    % solution for sink
    M_mu=RBModel.LF{ubi2}.ANq{n_mu}/5; % Z

    for kk=1:n_mu-1
        M_mu=M_mu+cond(kk)*RBModel.LF{ubi2}.ANq{kk}; % Stiffness
    end

    zN2=M_mu\RBModel.LF{ubi2}.FNq{1};
    zNh2 = RBModel.LF{ubi2}.V*zN2;

    Velecf1 = zNh1-zNh2;

    % Potential of the electrodes
    Velecf{rr}=Velecf1-Velecf1(ref); % Referenced to Cz
    % Removal of current injection electrodes, Cz and bad channels.
    Velecf{rr}([ubi1 ubi2],:)=[]; 

    % Load data of the first pattern:
    load([data '/Result' num2str(folder) '/Results/measurements/pattern_' num2str(ubi1) '.mat'], 'Data')
    u1 = Data.u(end-L+1:end,1);
    load([data '/Result' num2str(folder) '/Results/measurements/pattern_' num2str(ubi2) '.mat'], 'Data')
    u2 = Data.u(end-L+1:end,1);
    clear Data
    Vmeas1= u1-u2;

    % Remove bad channels
    % Note: Electrodes 6, 55 and 119 are always bad for all current injection
    % patterns in this dataset.
    Vmeas{rr}=Vmeas1-Vmeas1(ref);
    Vmeas{rr}([ubi1 ubi2],:)=[];
    
    Vmeas{rr} = abs(Vmeas{rr});
    u = Velecf{rr};

    plott = 1;
    % Figure of the electric potential at the electrodes
    if plott == 1
        figure; plot(1:size(u(:,:)),abs(u),1:size(u),Vmeas{rr});
        legend('ROM Forward Problem', 'Measurements');
        ylabel('Potential (module) [V]');
        disp(rr)
    end
    di = abs(Velecf{rr}-Vmeas{rr});
    [~,indx] = maxk(di,2);
    di(indx)=[];
    norms(rr,:) = norm(di);
    Vmeas2 = Vmeas{rr}; Vmeas2(indx)=[];
    Velecf2 = Velecf{rr}; Velecf2(indx)=[];
    rho(rr,:) = corr(Vmeas2,abs(Velecf2),'Type','Kendall');
    %pause(5)
end