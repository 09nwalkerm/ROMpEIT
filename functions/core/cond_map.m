function cond_map(p,t,f,estimates,electrodes)

figure()

t2 = t(t(:,5)==2,1:4);

elec_centers = zeros(size(estimates,1),3);
for ii = 1:length(electrodes)
    nodes = f(f(:,4)==electrodes(ii),1:3);
    nodes = unique(nodes);
    elec_centers(ii,:) = [mean(p(nodes,1)) mean(p(nodes,2)) mean(p(nodes,3))];
end

indx = unique(t2(:,1:4));
pp = p(indx,:);

PQ = dsearchn(pp,elec_centers);

F = scatteredInterpolant(pp(PQ,:),estimates(:,1));

VQ = F(p);

femeg_vis3d(p,t2,VQ,'p(:,1)>0')

colorbar()

end


























