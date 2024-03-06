function [FOM,RCOND]=femeg_ROM_RRBF_smart(FOM)
    RCOND = zeros(size(FOM.mu_train,1),1);
    for i=1:size(FOM.mu_train,1)
        NI = i;
        SRBF = SRBFunction(FOM,NI);
        RCOND(i) = rcond(SRBF);
        disp(RCOND(i))
%         if (RCOND == 0)
%             continue
%         elseif (RCOND < 1e-10)
%             SRBF = SRBFunction(FOM,i-1);
%             NI = i-1;
%             break
%         end
    end
    
    FOM.mu_train = FOM.mu_train(1:NI,:);
    FOM.betaa = FOM.betaa(1:NI);
    RRBF=SRBF\[log(FOM.betaa);zeros(length(FOM.active)+1,1)];
    FOM.RRBF=RRBF;
end

function SRBF = SRBFunction(FOM,NI)
    fi_i=@(x)exp(-x.^2) ;%x.^2 .* log(x);
    %fi_i=@(x)x.^2.*log(x);
    mu_train=FOM.mu_train(1:NI,:);
    nt=size(mu_train,1);
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
end