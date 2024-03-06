function [b,c,d,V,a] = femeg_vol_coord(p,t)
%FEMEG_VOL_COORD computes volume coordinates for a given mesh
%
%   [b,c,d,V,a] = FEMEG_VOL_COORD ( p, t) computes the volume coordinates a, b, 
%   c, d, and the volume V of each element of the mesh defined by p and t.
%
%   Inputs:
%      p: nodes defining the tetrahedral mesh (size: Np x 3).
%      t: elements defining the tetrahedral mesh (size: Nt x No, with No:
%         number of nodes for a given basis function, i.e. No=4 for first 
%         order basis functions, and No=10 for second order basis functions).
%
%   Outputs: 
%      a/b/c/d: Matrix with the corresponding coordinates for each element.
%               (size: Nt x 4, with e.g. a=[a1,a2,a3,a4], ai: ith volume
%               coordinate relative to a). "a" stands for the independent
%               term, "b" for the "x" term, "c" for the y term, and "d" for 
%               the z term. 
%      V: Element volumes (size: Nt x 1)
%
%
%  The notation is as used in the Appendix in:
%
%   + Beltrachini, L., "A finite element solution of the EEG forward problem for 
%      multipolar sources", Submitted

%   This file is part of the FEMEG toolbox.
%   Author: Leandro Beltrachini <BeltrachiniL at cardiff.ac.uk>

x1=p(t(:,1),1);x2=p(t(:,2),1);x3=p(t(:,3),1);x4=p(t(:,4),1);
y1=p(t(:,1),2);y2=p(t(:,2),2);y3=p(t(:,3),2);y4=p(t(:,4),2);
z1=p(t(:,1),3);z2=p(t(:,2),3);z3=p(t(:,3),3);z4=p(t(:,4),3);

b=[-y3.*z4+z3.*y4+y2.*z4-z2.*y4-y2.*z3+z2.*y3,y3.*z4-z3.*y4-y1.*z4+z1.*y4+y1.*z3-z1.*y3,-y2.*z4+z2.*y4+y1.*z4-z1.*y4-y1.*z2+z1.*y2,y2.*z3-z2.*y3-y1.*z3+z1.*y3+y1.*z2-z1.*y2];
c=[x3.*z4-x4.*z3-x2.*z4+x4.*z2+x2.*z3-x3.*z2,-x3.*z4+x4.*z3+x1.*z4-x4.*z1-x1.*z3+x3.*z1,x2.*z4-x4.*z2-x1.*z4+x4.*z1+x1.*z2-x2.*z1,-x2.*z3+x3.*z2+x1.*z3-x3.*z1-x1.*z2+x2.*z1];
d=[-x3.*y4+x4.*y3+x2.*y4-x4.*y2-x2.*y3+x3.*y2,x3.*y4-x4.*y3-x1.*y4+x4.*y1+x1.*y3-x3.*y1,-x2.*y4+x4.*y2+x1.*y4-x4.*y1-x1.*y2+x2.*y1,x2.*y3-x3.*y2-x1.*y3+x3.*y1+x1.*y2-x2.*y1];

if nargout==5
    a=[x2.*y3.*z4-x2.*y4.*z3-x3.*y2.*z4+x3.*y4.*z2+x4.*y2.*z3-x4.*y3.*z2,x1.*y4.*z3-x1.*y3.*z4+x3.*y1.*z4-x3.*y4.*z1-x4.*y1.*z3+x4.*y3.*z1,...
        x1.*y2.*z4-x1.*y4.*z2-x2.*y1.*z4+x2.*y4.*z1+x4.*y1.*z2-x4.*y2.*z1,x1.*y3.*z2-x1.*y2.*z3+x2.*y1.*z3-x2.*y3.*z1-x3.*y1.*z2+x3.*y2.*z1];
end

V=femeg_tetraquad_l(p,t(:,1:4));
    
end




function Q=femeg_tetraquad_l(p,t,N,F)

if nargin==2,N=1;end

% valid for each element (equally)
[q1,w1]=rquad(N,2);[q2,w2]=rquad(N,1);[q3,w3]=rquad(N,0);
[q1,q2,q3]=meshgrid(q1,q2,q3);q1=q1(:);q2=q2(:);q3=q3(:);
x=1-q1;y=(1-q2).*q1;z=q1.*q2.*q3;
w=abs(reshape(reshape(w2*w1',N^2,1)*w3',N^3,1));


% for each element
c=p(vec(t(:,1:4)'),:);
c(2:4:end,:)=c(2:4:end,:)-p(t(:,1),:);
c(3:4:end,:)=c(3:4:end,:)-p(t(:,1),:);
c(4:4:end,:)=c(4:4:end,:)-p(t(:,1),:);

% Wa=zeros(size(t,1),1);for k=1:size(t,1),Wa(k)=abs(det(c(2+(k-1)*4:4+(k-1)*4,:)));end

c2=reshape(c',12,size(t,1));
Wa=c2(4,:).*c2(8,:).*c2(12,:)-c2(4,:).*c2(11,:).*c2(9,:)-c2(7,:).*c2(5,:).*c2(12,:)...
    +c2(7,:).*c2(11,:).*c2(6,:)+c2(10,:).*c2(5,:).*c2(9,:)-c2(10,:).*c2(8,:).*c2(6,:);

if nargin==2
    Q=Wa*sum(w);
else
    A=[ones(N^3,1) x y z];
    
    X=(A*reshape(c(:,1),4,size(t,1)))';
    Y=(A*reshape(c(:,2),4,size(t,1)))';
    Z=(A*reshape(c(:,3),4,size(t,1)))';
    
    Q=Wa.*(w'*feval(F,X,Y,Z)');
end
Q=Q';

end


function [x,w]=rquad(N,k)

k1=k+1; k2=k+2; 
n=1:N;  nnk=2*n+k;
A=[k/k2 repmat(k^2,1,N)./(nnk.*(nnk+2))];
n=2:N; nnk=nnk(n);
B1=4*k1/(k2*k2*(k+3)); nk=n+k; nnk2=nnk.*nnk;
B=4*(n.*nk).^2./(nnk2.*nnk2-nnk2);
ab=[A' [(2^k1)/k1; B1; B']]; s=sqrt(ab(2:N,2));
[V,X]=eig(diag(ab(1:N,1),0)+diag(s,-1)+diag(s,1));
[X,I]=sort(diag(X));    

% Grid points
x=(X+1)/2; 

% Quadrature weights
w=(1/2)^(k1)*ab(1,2)*V(1,I)'.^2;

end


function v=vec(a)

v=a(:);

end