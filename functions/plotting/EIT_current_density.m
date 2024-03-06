function J=EIT_current_density(p,t,D,u)

[b,c,d,Vtot]=femeg_vol_coord(p,t(:,1:4));

J=[dot(repmat(D(:,1),1,4).*b+repmat(D(:,2),1,4).*c+repmat(D(:,3),1,4).*d,u(t(:,1:4)),2)./Vtot,...
   dot(repmat(D(:,2),1,4).*b+repmat(D(:,4),1,4).*c+repmat(D(:,5),1,4).*d,u(t(:,1:4)),2)./Vtot,...
   dot(repmat(D(:,3),1,4).*b+repmat(D(:,5),1,4).*c+repmat(D(:,6),1,4).*d,u(t(:,1:4)),2)./Vtot];

J=-J/6;

end