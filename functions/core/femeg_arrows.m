function [CBM,Dn] = femeg_arrows(rank)

load head_model_theta

tissue = 2; %skull
t2 = t(t(:,5)==tissue);
nt = length(t2(:,1));
oo = zeros(length(t2(:,1)),1);
% for ii=1:length(theta(:,1))
%     if theta(ii,1) 
%         theta(:,1) = theta(:,1) + pi/2;
%theta(:,:) = -1.*theta(:,:);

%CBM = [cos(theta(:,1)) cos(theta(:,2)).*sin(theta(:,1)) sin(theta(:,2)).*sin(theta(:,1)); -sin(theta(:,1)) cos(theta(:,2)).*cos(theta(:,1)) sin(theta(:,2)).*cos(theta(:,1)); oo -sin(theta(:,2)) cos(theta(:,2))];
%CBM = [cos(theta(:,1)) -sin(theta(:,1)) oo cos(theta(:,2)).*sin(theta(:,1)) cos(theta(:,2)).*cos(theta(:,1)) -sin(theta(:,2)) sin(theta(:,2)).*sin(theta(:,1)) sin(theta(:,2)).*cos(theta(:,1)) cos(theta(:,2))];
%CBM = [sin(theta(:,2)).*cos(theta(:,1)) sin(theta(:,2)).*sin(theta(:,1)) cos(theta(:,2)) cos(theta(:,2)).*cos(theta(:,1)) cos(theta(:,2)).*sin(theta(:,1)) -sin(theta(:,2)) -sin(theta(:,1)) cos(theta(:,1)) oo];
%CBM = [cos(theta(:,2)) cos(theta(:,1)).*sin(theta(:,2)) sin(theta(:,2)).*sin(theta(:,1)) -sin(theta(:,2)) cos(theta(:,2)).*cos(theta(:,1)) sin(theta(:,1)).*cos(theta(:,2)) oo -sin(theta(:,1)) cos(theta(:,1))];
%CBM = [cos(theta(:,1)) -sin(theta(:,1)) oo sin(theta(:,1)) cos(theta(:,1)) oo oo oo oo+1];
%CBM = [sin(theta(:,2)).*cos(theta(:,1)) cos(theta(:,2)).*cos(theta(:,1)) -sin(theta(:,1)) sin(theta(:,2)).*sin(theta(:,1)) cos(theta(:,2)).*sin(theta(:,1)) cos(theta(:,1)) cos(theta(:,2)) -sin(theta(:,2)) oo];
%CBM = [sin(theta(:,1)).*cos(theta(:,2)) cos(theta(:,1)).*cos(theta(:,2)) -sin(theta(:,2)) sin(theta(:,1)).*sin(theta(:,2)) cos(theta(:,1)).*sin(theta(:,2)) cos(theta(:,2)) cos(theta(:,1)) -sin(theta(:,1)) oo];
%CBM = [cos(theta(:,1)).*cos(theta(:,2)) sin(theta(:,1)).*cos(theta(:,2)) -sin(theta(:,2)) -sin(theta(:,1)) cos(theta(:,1)) oo cos(theta(:,1)).*sin(theta(:,2)) sin(theta(:,1)).*sin(theta(:,2)) cos(theta(:,2))];
CBM = [cos(theta(:,1)).*cos(theta(:,2)) sin(theta(:,1)).*cos(theta(:,2)) -sin(theta(:,2)) -sin(theta(:,1)) cos(theta(:,1)) oo cos(theta(:,1)).*sin(theta(:,2)) sin(theta(:,1)).*sin(theta(:,2)) cos(theta(:,2))];
%CBM = [cos(theta(:,1)).*cos(theta(:,2)) -sin(theta(:,1)) cos(theta(:,1)).*sin(theta(:,2)) sin(theta(:,1)).*cos(theta(:,2)) cos(theta(:,1)) sin(theta(:,1)).*sin(theta(:,2)) -sin(theta(:,2)) oo cos(theta(:,2))]; %Transpose of above

Dn = zeros(nt,6);Dn(:,6)=0.01;Dn(:,[1,4])=0.005;
%Dn(1:10,:)
Dn = change_basis(Dn,CBM,"all");

%[t3,index] = unique(t2(:,1));
%Dn = Dn(index,:);

% X = p(t3(:,1),1);
% Y = p(t3(:,1),2);
% Z = p(t3(:,1),3);
% % psi = Dn(:,3);
% % theta2 = Dn(:,5);
% % R = Dn(:,6);
% % U = R.*sin(theta2).*cos(psi);
% % V = R.*sin(theta2).*sin(psi);
% % W = R.*cos(theta2);
% U = Dn(:,3);
% V = Dn(:,5);
% W = Dn(:,6);
% 
% quiver3(X,Y,Z,U,V,W)
scale = 0.01;
figure()

% USE THIS ONE
if rank == 2
    for ii=1:30:30000
        D = [Dn(ii,1) Dn(ii,2) Dn(ii,3);Dn(ii,2) Dn(ii,4) Dn(ii,5); Dn(ii,3) Dn(ii,5) Dn(ii,6)];
        %disp(D)
        [vec,val] = eig(D);
        %disp(vec)
        ee = diag(val);
        [~,ind] = sort(diag(val),'descend');
        X = centroids(ii,1);
        Y = centroids(ii,2);
        Z = centroids(ii,3);
        U = vec(1,ind(1));%Dn(ii,1);%
        V = vec(2,ind(1));%Dn(ii,2);%
        W = vec(3,ind(1));%Dn(ii,3);%
        quiver3(X,Y,Z,U,V,W,ee(ind(1)))
        hold on
%         U = vec(1,ind(2));
%         V = vec(2,ind(2));
%         W = vec(3,ind(2));
%         quiver3(X,Y,Z,U,V,W,scale)
%         U = vec(1,ind(3));
%         V = vec(2,ind(3));
%         W = vec(3,ind(3));
%         quiver3(X,Y,Z,U,V,W,scale)
    end
elseif rank==1
    for ii=1:30:30000
        X = centroids(ii,1);
        Y = centroids(ii,2);
        Z = centroids(ii,3);
        U = Dn(ii,3);
        V = Dn(ii,5);
        W = Dn(ii,6);
        quiver3(X,Y,Z,U,V,W,scale)
        hold on
    end
elseif rank == 3
    for ii=1:30:30000
        X = centroids(ii,1);
        Y = centroids(ii,2);
        Z = centroids(ii,3);
        U = CBM(ii,1);%Dn(ii,1);%
        V = CBM(ii,2);%Dn(ii,2);%
        W = CBM(ii,3);%Dn(ii,3);%
        quiver3(X,Y,Z,U,V,W,scale)
        hold on
        U = CBM(ii,4);
        V = CBM(ii,5);
        W = CBM(ii,6);
        quiver3(X,Y,Z,U,V,W,scale)
        U = CBM(ii,7);
        V = CBM(ii,8);
        W = CBM(ii,9);
        quiver3(X,Y,Z,U,V,W,scale)
    end
end
% hold on
% 
% U = Dn(:,1);
% V = Dn(:,2);
% W = Dn(:,3);
% 
% quiver3(X,Y,Z,U,V,W)
% 
% U = Dn(:,2);
% V = Dn(:,4);
% W = Dn(:,5);
% 
% quiver3(X,Y,Z,U,V,W)

end