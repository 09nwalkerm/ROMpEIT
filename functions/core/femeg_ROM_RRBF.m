function FOM=femeg_ROM_RRBF(FOM)

fi_i=@(x)exp(-x.^2) ;%x.^2 .* log(x);
%fi_i=@(x)x.^2.*log(x);
nt=size(FOM.mu_train,1);
mu_train=FOM.mu_train;
M=zeros(nt,nt);
for ii=1:nt
    for jj=1:nt
        M(ii,jj)=fi_i(norm(mu_train(ii,:)-mu_train(jj,:)));
    end
end

% It can be seen that P=mu_train';
osm=ones(nt,1);
%SRBF=[M,mu_train,osm;mu_train',zeros(length(FOM.active),length(FOM.active)+1);osm',zeros(1,length(FOM.active)+1)];
SRBF=[M,mu_train,osm;mu_train',zeros(length(FOM.active),length(FOM.active)+1);osm',zeros(1,length(FOM.active)+1)];
RRBF=SRBF\[log(FOM.betaa);zeros(length(FOM.active)+1,1)];
FOM.RRBF=RRBF;