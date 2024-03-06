function [Dn] = change_basis(Dn,CBM,direction)

%% Change basis function

% This function changes the basis of the conductivity tensor to a spherical
% one with tangential and radial components. It uses the change of basis
% matrix (CMB) and the conductivity tensor matrix (Dn)

% Conductivity tensor in the form: [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz]
   %                               [D1 ,D2 ,D3 ,D4 ,D5 ,D6 ]

%% function

%CBM_inv = inv(CBM);
D = Dn;
% vectorise
%CBM = CBM(:)';
%fprintf('Olo')
%CBM_inv = vec(CBM_inv)';

% Calc CMB*Dn
if direction == "tan"
    % Calc CMB*Dn
    %DNN = [CBM(:,1).*Dn(:,1),CBM(:,2).*Dn(:,1),0,0,0,CBM(:,6).*Dn(:,5),0,0,0];
    % Clc DNN*CMB'
    % Dn = [DNN(:,1).*CMB_inv(:,1),DNN(:,2).*CMB_inv(:,1),0,DNN(:,1).*CMB_inv(:,4),DNN(:,2).*CMB_inv(:,4),0,0,0,DNN(:,6).*CMB_inv(:,8)];
    %Dn(:,1) = DNN(:,1).*CBM_inv(:,1);
    %Dn(:,2) = DNN(:,2).*CBM_inv(:,1);
    %Dn(:,4) = DNN(:,1).*CBM_inv(:,4);
    %Dn(:,5) = DNN(:,2).*CBM_inv(:,4);
    %Dn(:,9) = DNN(:,6).*CBM_inv(:,8);
    %Dn(:,1) = DNN(:,1);
    %Dn(:,2) = DNN(:,2);
    %Dn(:,6) = DNN(:,6);
    Dn(:,1) = CBM(:,1).^2.*D(:,1) + CBM(:,4).^2.*D(:,4);
    Dn(:,2) = CBM(:,1).*D(:,1).*CBM(:,2) + CBM(:,4).*D(:,4).*CBM(:,5);
    Dn(:,3) = CBM(:,4).*D(:,4).*CBM(:,6) + CBM(:,1).*D(:,1).*CBM(:,3);
    Dn(:,4) = CBM(:,2).^2.*D(:,1) + CBM(:,5).^2.*D(:,4);
    Dn(:,5) = CBM(:,5).*D(:,4).*CBM(:,6) + CBM(:,2).*D(:,1).*CBM(:,3);
    Dn(:,6) = CBM(:,6).^2.*D(:,4) + CBM(:,3).^2.*D(:,1);
    %fprintf('Olotan')
    
elseif direction == "rad"
    % Calc CMB*Dn
    %DNN = [0,0,0,0,0,0,CBM(:,7).*Dn(:,9),CBM(:,8).*Dn(:,9),0];
    % Clc DNN*CMB'
    % Dn = [DNN(:,7).*CMB_inv(:,3),DNN(:,8).*CMB_inv(:,3),0,0,0,0,0,0,0];
    %Dn(:,1) = DNN(:,7).*CBM_inv(:,3);
    %Dn(:,2) = DNN(:,8).*CBM_inv(:,3);
    %Dn(:,7) = DNN(:,7);
    %Dn(:,8) = DNN(:,8);
    Dn(:,1) = CBM(:,7).^2.*D(:,6);
    Dn(:,2) = CBM(:,7).*D(:,6).*CBM(:,8);
    Dn(:,3) = CBM(:,7).*D(:,6).*CBM(:,9);
    Dn(:,4) = CBM(:,8).^2.*D(:,6);
    Dn(:,5) = CBM(:,8).*D(:,6).*CBM(:,9);
    Dn(:,6) = CBM(:,9).^2.*D(:,6);
    %fprintf('Olorad')
    
elseif direction == "all"
    Dn(:,1) = CBM(:,1).^2.*D(:,1) + CBM(:,4).^2.*D(:,4) + CBM(:,7).^2.*D(:,6);
    Dn(:,2) = CBM(:,1).*D(:,1).*CBM(:,2) + CBM(:,4).*D(:,4).*CBM(:,5) + CBM(:,7).*D(:,6).*CBM(:,8);
    Dn(:,3) = CBM(:,4).*D(:,4).*CBM(:,6) + CBM(:,7).*D(:,6).*CBM(:,9) + CBM(:,1).*D(:,1).*CBM(:,3);
    Dn(:,4) = CBM(:,2).^2.*D(:,1) + CBM(:,5).^2.*D(:,4) + CBM(:,8).^2.*D(:,6);
    Dn(:,5) = CBM(:,5).*D(:,4).*CBM(:,6) + CBM(:,8).*D(:,6).*CBM(:,9) + CBM(:,2).*D(:,1).*CBM(:,3);
    Dn(:,6) = CBM(:,6).^2.*D(:,4) + CBM(:,9).^2.*D(:,6) + CBM(:,3).^2.*D(:,1);
    fprintf('Oloall')
    
elseif direction == "vector"
    Dn(:,3) = CBM(:,7).*Dn(:,6);
    Dn(:,5) = CBM(:,8).*Dn(:,6);
    Dn(:,6) = CBM(:,9).*Dn(:,6);
    
end
    
end
