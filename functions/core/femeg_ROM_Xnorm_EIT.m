function [X]=femeg_ROM_Xnorm_EIT(p,t,L,c)
% This function computes the scalar product matrix defined in H^1_\ast for
% the EEG forward problem, as given in equation XX in RR.  

np=size(p,1);
[~,~,~,V,~]=femeg_vol_coord(p,t(:,1:4));


X1=femeg_ROM_Xnorm_EEG(p,t);


Mx=accumarray(t(:,1),V/4,[np,1],[],[])+accumarray(t(:,2),V/4,[np,1],[],[])...
  +accumarray(t(:,3),V/4,[np,1],[],[])+accumarray(t(:,4),V/4,[np,1],[],[]);
[i,j]=find(X1~=0);
X2=sparse(i,j,Mx(i)+Mx(j),np,np);

B1=X1-c*X2+2*c^2*(X1~=0);

X=blkdiag(B1,speye(L));
