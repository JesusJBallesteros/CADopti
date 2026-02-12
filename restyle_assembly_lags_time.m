function[As_restyled,As_restyled_index]=restyle_assembly_lags_time(A,A_index)
%% REORDER LAGS AND SHIFT ASSEMBLY'S OCCURRENCE
% lags and elements are reordered so that all lags are expressed with respect 
% of the first firing unit. Assembly activation times is shifted
% accordingly to match the firing of the first assembly unit.
%
%
    shiftedbins=actbins+minlag;
    validbins=(shiftedbins>=1) & (shiftedbins<=numel(A{i}.Time));
    As_restyled{i}.Time(shiftedbins(validbins))=double(A{i}.Time(actbins(validbins)));
%  for information please contact eleonora.russo@zi-mannheim.de


As_restyled=cell(size(A));
As_restyled_index=A_index;


numAss=length(As_restyled);
for i=1:numAss
    llag=[0 A{i}.lag];
    [a b]=sort(llag);
    minlag=min(a);
    As_restyled{i}.elements=A{i}.elements(b);
    As_restyled{i}.lag=a-minlag; 

    aus=double(A{i}.Time);
    actbins=find(~~aus);

    As_restyled{i}.Time=zeros(size(A{i}.Time));
    As_restyled{i}.Time(actbins+minlag)= double(A{i}.Time(actbins));

    
    As_restyled{i}.pr=A{i}.pr(end);
    As_restyled{i}.Noccurrences=A{i}.Noccurrences(end);
    As_restyled{i}.bin=A{i}.bin;
    
end




end
