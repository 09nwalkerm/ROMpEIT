function [S]=prepare_EIT_model_cem(p,t,L,np)

list_layers=unique(t(:,5));
N_layers=length(list_layers); % get number of layers

% Compute the compartamental stiffness matrices
S=struct([]);
for kk=1:N_layers
    Ind=t(:,5)==list_layers(kk);
    DD=zeros(size(t,1),6);DD(Ind,[1,4,6])=1;
    S(kk).stiff=sparse([femeg_stiffness(p,t,DD),zeros(np,L);zeros(L,np+L)]);
end