function [assemD]=TestPair_ref(ensemble,spikeTrain2,n2,maxlag,Dc,reference_lag)
% this function tests if the two spike trains have repetitive patterns occurring more
% frequently than chance.

% ensemble := structure with the previously formed assembly and its spike train
% spikeTrain2 := spike train of the new unit to be tested for significance (candidate to be a new assembly member)
% n2 := new unit tested
% maxlag := maximum lag to be tested
% Dc := size (in bins) of chunks in which I divide the spike trains to compute the variance (to reduce non stationarity effects on variance estimation)
% reference_lag := lag of reference; if zero or negative reference lag=-l
%
%  � 2020 Russo
%  for information please contact eleonora.russo@zi-mannheim.de

couple=[double(ensemble.Time)-min(double(ensemble.Time(:))) double(spikeTrain2)-min(double(spikeTrain2(:)))];
nu=2;
ntp=size(couple,1);

maxrate=max(couple(:));
if maxrate==0
    ExpAB=0;
else
    marg_all=marginal_counts_by_level(couple,maxrate);
    ExpAB=sum((marg_all(1,:).*marg_all(2,:))/ntp);
end

if reference_lag>maxlag
    ctAB=nan(1,reference_lag+1);
    ctAB_=nan(1,reference_lag+1);
    for l=0:reference_lag
        trAB=[couple(1:end-reference_lag,1), couple(l+1:end-reference_lag+l,2)];
        trBA=[couple(l+1:end-reference_lag+l,1), couple(1:end-reference_lag,2)];
        ctAB(l+1)=sum(min(trAB'));
        ctAB_(l+1)=sum(min(trBA'));
    end
    Hab_l=[ctAB_(end:-1:2), ctAB];
    [a,b]=max(Hab_l(reference_lag-maxlag+1:end-reference_lag+maxlag));
    lags=-reference_lag:reference_lag;
    l_=lags(b+reference_lag-maxlag);
    Hab=a;
    if l_<0, l_ref=l_+reference_lag; else, l_ref=l_-reference_lag; end
    idx_ref=find(lags==l_ref,1);
    if isempty(idx_ref)
        l_ref=max(min(l_ref,lags(end)),lags(1));
        idx_ref=find(lags==l_ref,1);
    end
    Hab_ref=Hab_l(idx_ref);
else
    ctAB=nan(1,maxlag+1);
    ctAB_=nan(1,maxlag+1);
    for l=0:maxlag
        trAB=[couple(1:end-maxlag,1), couple(l+1:end-maxlag+l,2)];
        trBA=[couple(l+1:end-maxlag+l,1), couple(1:end-maxlag,2)];
        ctAB(l+1)=sum(min(trAB'));
        ctAB_(l+1)=sum(min(trBA'));
    end

    if reference_lag<=0
        aus=[ctAB; ctAB_];
        [a,b]=max(aus(:));
        [I,J]=ind2sub(size(aus),b);
        l_=(I==1)*(J-1)-(I==2)*(J-1);
        if l_==0
            Hab=ctAB(1);
            if numel(ctAB_)>=3, Hab_ref=ctAB_(3); else, Hab_ref=ctAB_(1); end
        else
            Hab=ctAB(J); Hab_ref=ctAB_(J);
        end
    else
        Hab_l=[ctAB_(end:-1:2), ctAB];
        [a,b]=max(Hab_l(:));
        lags=-maxlag:maxlag;
        l_=lags(b);
        Hab=a;
        if l_<0, l_ref=l_+reference_lag; else, l_ref=l_-reference_lag; end
        idx_ref=find(lags==l_ref,1);
        if isempty(idx_ref)
            l_ref=max(min(l_ref,lags(end)),lags(1));
            idx_ref=find(lags==l_ref,1);
        end
        Hab_ref=Hab_l(idx_ref);
    end
end
TPrMTot=[0 Hab; Hab_ref 0];

if  a==0  || ExpAB<=5  || ExpAB>=(min(sum(couple))-5)
    assemD.elements=uint32([double(ensemble.elements) n2]);
    assemD.lag=int16([double(ensemble.lag), 99]);
    assemD.pr=[double(ensemble.pr) 1];
    assemD.Time=[];
    assemD.Noccurrences=uint32([double(ensemble.Noccurrences) 0]);
else
    len=size(couple,1);
    if l_==0
        Time=min(couple(:,1),couple(:,2));
    elseif l_>0
        Time=zeros(len,1);
        Time(1:len-l_)=min(couple(1:end-l_,1),couple(l_+1:end,2));
    else
        Time=zeros(len,1);
        Time(-l_+1:end)=min(couple(-l_+1:end,1),couple(1:end+l_,2));
    end

    nch=ceil((ntp-maxlag)/Dc);
    Dc=floor((ntp-maxlag)/nch);

    n=ntp-maxlag;
    varXtot=zeros(2);

    for iii=1:nch
        st=(1+Dc*(iii-1));
        if iii<nch, en=Dc*iii; else, en=(ntp-maxlag); end

        if l_==0
            couple_t=couple(st:en,:);
        elseif l_>0
            couple_t=[couple(st:en,1), couple(st+l_:en+l_,2)];
        else
            couple_t=[couple(st-l_:en-l_,1), couple(st:en,2)];
        end

        ch_n=size(couple_t,1);
        maxrate_t=max(couple_t(:));
        marg_pr=marginal_counts_by_level(couple_t,maxrate_t);

        varT=zeros(nu);
        covX=zeros(nu);
        for i=1:maxrate_t
            Mx_i=marg_pr(:,i)*ones(1,2);
            cov_i=(Mx_i.*Mx_i'./ch_n).*(ch_n-Mx_i).*(ch_n-Mx_i')./(ch_n*(ch_n-1));
            varT=varT+cov_i;
            covX=covX+cov_i/(ch_n-1);

            for j=(i+1):maxrate_t
                Mx_j=marg_pr(:,j)*ones(1,2);
                cov_ij=2*(Mx_j.*Mx_j'./ch_n).*(ch_n-Mx_i).*(ch_n-Mx_i')./(ch_n*(ch_n-1));
                varT=varT+cov_ij;
                covX=covX+cov_ij/(ch_n-1);
            end
        end

        varX=varT+varT'-covX-covX';
        varXtot=varXtot+varX;
    end

    X=TPrMTot-TPrMTot';
    if abs(X(1,2))>0
        X=abs(TPrMTot-TPrMTot')-0.5;
    end

    if varXtot(1,2)==0
        prF=1;
    else
        F=X.^2./varXtot;
        prF=fcdf(F(1,2),1,n,'upper');
    end

    assemD.elements=uint32([double(ensemble.elements) n2]);
    assemD.lag=int16([double(ensemble.lag), l_]);
    assemD.pr=[double(ensemble.pr) prF];
    assemD.Time=uint16(max(0,Time));
    assemD.Noccurrences=uint32([double(ensemble.Noccurrences) sum(Time)]);
end

end

function marg_pr=marginal_counts_by_level(couple,maxrate)
if maxrate<=0
    marg_pr=zeros(2,0);
    return
end
v1=accumarray(couple(:,1)+1,1,[maxrate+1,1]);
v2=accumarray(couple(:,2)+1,1,[maxrate+1,1]);
c1=flipud(cumsum(flipud(v1)));
c2=flipud(cumsum(flipud(v2)));
marg_pr=[c1(2:end)'; c2(2:end)'];
end