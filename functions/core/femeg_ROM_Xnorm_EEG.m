function X=femeg_ROM_Xnorm_EEG(p,t)
% This function computes the scalar product matrix defined in H^1_\ast for
% the EEG forward problem, as given in equation 3.26??? in Somersalo et al. (1992)  

np=size(p,1);
[~,~,~,V,~]=femeg_vol_coord(p,t(:,1:4));

% Compute first term
Mv=(ones(4)+eye(4))/120;
Me=6*V*Mv(:)';
row_id=repmat(t(:,1:4),1,4);
col_id=0*row_id;cont=0;
for kk=1:4:4*4,cont=cont+1;col_id(:,kk:kk+4-1)=repmat(t(:,cont),1,4);end
X=accumarray([row_id(:,1),col_id(:,1)],Me(:,1),[np,np],[],[],true);
for k=2:size(Mv,1)^2
    X=X+accumarray([row_id(:,k),col_id(:,k)],Me(:,k),[np,np],[],[],true);
end

% Sum up the second term
X=X+femeg_stiffness(p,t,repmat([1,0,0,1,0,1],size(t,1),1)); %Actually the first term ||del u ||^2_L_2
