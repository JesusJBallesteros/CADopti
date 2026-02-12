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


couple=[ensemble.Time-min(ensemble.Time(:)) spikeTrain2-min(spikeTrain2(:))]; % spike train pair I am going to test
nu=2;
ntp=size(couple,1);  %%%% trial length

maxrate=max(couple(:));
ExpABi=zeros(1,maxrate);
for i=1:maxrate
    level_hits=(couple>=i);
    ExpABi(i)=prod(sum(level_hits,1))/ntp;
end

if reference_lag>maxlag
    % % % decide which is the lag with most coincidences (l_:=best lag)
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
    if l_<0
        l_ref=l_+reference_lag;
        Hab_ref=Hab_l(find(lags==l_ref));
    else
        l_ref=l_-reference_lag;
        Hab_ref=Hab_l(find(lags==l_ref));
    end
else
    % % % decide which is the lag with most coincidences (l_:=best lag)
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
        [I,J] = ind2sub(size(aus),b);
        l_=(I==1)*(J-1)-(I==2)*(J-1);  %% I select l_
        if l_==0
            Hab=ctAB(1); Hab_ref=ctAB_(3);
        else
            Hab=ctAB(J); Hab_ref=ctAB_(J);
        end
    else
        Hab_l=[ctAB_(end:-1:2), ctAB];
        [a,b]=max(Hab_l(:));
        lags=-maxlag:maxlag;
        l_=lags(b);
        Hab=a;
        if l_<0
            l_ref=l_+reference_lag;
            Hab_ref=Hab_l(find(lags==l_ref));
        else
            l_ref=l_-reference_lag;
            Hab_ref=Hab_l(find(lags==l_ref));
        end

    end

end
TPrMTot=[0 Hab; Hab_ref 0]; % matrix with #AB and #BA


%% HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH


ExpAB=sum(ExpABi);

if  a==0  || ExpAB<=5  || ExpAB>=(min(sum(couple))-5)    %case of no coincidences or limit for the F asinptotical distribution (too few coincidences)
    assemD.elements=[ensemble.elements n2];
    assemD.lag=[ensemble.lag, 99];
    assemD.pr=[ensemble.pr 1];  % setting pr=1 the tested pair will be discarded as assembly
    assemD.Time=[];
    assemD.Noccurrences=[ensemble.Noccurrences 0];
else

    % % % construct the activation time series for the couple
    len=size(couple,1);        %%%% trial length
    Time=uint8(zeros(len,1));  %%%% activation vector

    if l_==0
        for i=1:maxrate
            level_hits=(couple>=i);
            Time=Time+level_hits(:,1).*level_hits(:,2);
        end
    elseif l_>0
        for i=1:maxrate
            level_hits=(couple>=i);
            Time(1:len-l_)=Time(1:len-l_)+level_hits(1:end-l_,1).*level_hits(l_+1:end,2);
        end
    else
        for i=1:maxrate
            level_hits=(couple>=i);
            Time(-l_+1:end)=Time(-l_+1:end)+level_hits(-l_+1:end,1).*level_hits(1:end+l_,2);
        end
    end


    %% --------------------------------------------------------------------%
    % % % HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    % I cut the spike train in stationary segments
    %%%%%%%%%%%%%%%%%%%%% chunking  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nch=ceil((ntp-maxlag)/Dc);
    Dc=floor((ntp-maxlag)/nch); %% new chunk size, this is to have all chunks of rougly the same size

    n=ntp-maxlag;
    varXtot=zeros(2);

    for iii=1:nch
        st=(1+Dc*(iii-1));
        if iii<nch
            en=Dc*iii;
        else
            en=(ntp-maxlag);
        end

        if l_==0
            couple_t=couple(st:en,:);
        elseif l_>0
            couple_t=[couple(st:en,1), couple(st+l_:en+l_,2)];
        else
            couple_t=[couple(st-l_:en-l_,1), couple(st:en,2)];
        end

        ch_n=size(couple_t,1);
        maxrate_t=max(couple_t(:));

        marg_pr=zeros(2,maxrate_t);
        for i=1:maxrate_t
            level_hits=(couple_t>=i);
            marg_pr(:,i)=[sum(level_hits(:,1)); sum(level_hits(:,2))];
        end

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
        X=abs(TPrMTot-TPrMTot')-0.5;   %Yates correction
    end

    if varXtot(1,2)==0  % if variance is zero
        prF=1;
    else
        F=X.^2./varXtot;
        prF=fcdf(F(1,2),1,n,'upper');
    end



    %%
    %All information about the assembly and test are returned

    assemD.elements=[ensemble.elements n2];
    assemD.lag=[ensemble.lag, l_];
    assemD.pr=[ensemble.pr prF];
    assemD.Time=Time;
    assemD.Noccurrences=[ensemble.Noccurrences sum(Time)];


end

end
