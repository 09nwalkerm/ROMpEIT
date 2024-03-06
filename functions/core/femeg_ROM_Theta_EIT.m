function [ theta_a, theta_f ] = femeg_ROM_Theta_EIT( mu, ROM )

theta_a=[mu(1:end-1),1/mu(end)];
theta_f=1;

%if ROM.P == ROM.active(end)
%    theta_a=[mu(1:end-1),1/mu(end)];
%    theta_f=1;
%else
%    theta_a=mu;
%    theta_f=1;
%end
