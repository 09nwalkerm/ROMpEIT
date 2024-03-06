function [varargout]=femeg_stiffness_cem(p,t,f,z,varargin)
%FEMEG_STIFFNESS_CEM computes the stiffness matrix of the EEG-FP for the
%                    complete electrode model
%
%   FEMEG_STIFFNESS_CEM ( p, t, f, z) computes the stiffness matrix of the EEG
%   forward problem considering the CEM given the mesh with nodes p, elements t, 
%   and external faces f. This is done as described in:  
%
%   + Beltrachini, L., "MY CEM", Submitted
%
%
%   Inputs:
%      p: nodes defining the tetrahedral mesh (size: Np x 3).
%      t: elements defining the tetrahedral mesh (size: Nt x No, with No:
%         number of nodes for a given basis function, i.e. No=4 for first 
%         order basis functions, and No=10 for second order basis functions).
%      f
%      z:
%      S (optional): Stiffness matrix as for the PEM as computed by
%        femeg_stiffness. If provided, it will also compute the final 
%        CEM stiffness matrix (size: Np+L x Np+L).
%
%   Outputs:
%      If S was provided, it will return:
%      
%      Scem= Ensembled CEM stiffness matrix (size: (Np+L) x (Np+L)).
%      A: complemetary matrix to the PEM stiffness matrix (size: Np x Np).
%      B: (size: Np x L)
%      C: (size: L x L)
%
%      If S is not provided, it will return A, B, and C only. 
%
% The function femeg_som can be used for defining a second order mesh as a
% function of a first order mesh.
%
% See also FEMEG_SOM

%   This file is part of the FEMEG toolbox.
%   Author: Leandro Beltrachini <BeltrachiniL at cardiff.ac.uk>


% Preliminaries
np=size(p,1); % number of nodes
st=size(t,2);if st<=5,order=1;nel=3;else order=2;nel=6;end % get FEM order
L=length(unique(f(:,end)))-1; % number of electrodes
if length(z)==1,z=z*ones(L,1);end
A=spalloc(np,np,50*np);
B=spalloc(np,L,50*L);
C=spalloc(L,L,L);

% Compute external triangles area
areaa=cross(p(f(:,2),:)-p(f(:,1),:),p(f(:,1),:)-p(f(:,3),:),2);
areaa=sqrt(sum(areaa.^2,2))/2;

% Iterate for each electrode 
for kk=1:L
    
    Ind=f(:,end)==kk;
    
    %%% Compute matrix A
    if order==1
        A_2e=(eye(3)+ones(3))/24;
    elseif order==2
        A_2e=[6,-1,-1,0,0,-4;-1,6,-1,-4,0,0;-1,-1,6,0,-4,0;...
              0,-4,0,32,16,16;0,0,-4,16,32,16;-4,0,0,16,16,32]/360;
    end
    A_2e=2/z(kk)*areaa(Ind)*A_2e(:)';
    
    % Ensamble matrix
    tt=f(Ind,1:nel);
    row_id=repmat(tt(:,1:nel),1,nel);
    col_id=kron(tt(:,1:nel),ones(1,nel));
    
    A_2n=accumarray([row_id(:,1),col_id(:,1)],A_2e(:,1),[np,np],[],[],true);
    for jj=2:size(A_2e,2)
        A_2n=A_2n+accumarray([row_id(:,jj),col_id(:,jj)],A_2e(:,jj),[np,np],[],[],true);
    end
    A=A+A_2n;
     
    %%% Compute matrix B
    if order==1
        B_e=1/z(kk)*areaa(Ind)*[1,1,1]/3;
    elseif order==2
        B_e=1/z(kk)*areaa(Ind)*[0,0,0,1,1,1]/3;
    end
    B=B+accumarray([tt(:),kk*ones(nel*size(tt,1),1)],B_e(:),[np,L],[],[],true);

    %%% Compute matrix C
    C=C+sparse(kk,kk,sum(areaa(Ind))/z(kk),L,L);
end

%%% If the PEM stiffness matrix was provided, ensemble the full CEM matrix
if nargin==5
    varargout{1}=[varargin{1}+A,-B;-B',C];
    varargout{2}=A;
    varargout{3}=B;
    varargout{4}=C;
else
    varargout{1}=A;
    varargout{2}=B;
    varargout{3}=C;
end
