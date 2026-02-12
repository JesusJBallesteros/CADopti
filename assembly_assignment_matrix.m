function [Amatrix,Binvector,Unit_order,As_order]=assembly_assignment_matrix(As_across_bins,nneu,BinSizes,display)
% Display the Assembly assignment matrix.
% USAGE: [Amatrix,Binvector,Unit_order]=assembly_assignment_matrix(As_across_bins,nneu,BinSizes,display)
%
% ARGUMENTS:
%   As_across_bins  - structure containing assembly information
%   nneu            - # of units
%   BinSizes          - list of investigated bin sizes
%   display         - Possible display options:
%           1) display='raw'  ->  no rearrangement made, assemblies and units displayed as 
%                                 from the detection algorithm
%           2) display='ordunit'  ->  Assembly order is left unchanged, units rearranged 
%                                     in order to separate assembly units from units not taking part to any assembly.
%           3) display='clustered'  -> Assembly and unit order changed in order to
%                                      display close to one other assemblies with similar elements
%
% RETURNS:
%       Amatrix     - ASSEMBLY ASSIGNMENT MATRIX: each column is an assembly 
%       each raw a unit. If unit i takes part in assembly j, Amatrix(i,j) displays
%       the firing delay (lag) of i with respect to the first assembly unit.
%       Binvector   - bin size. Binvector(j) is the bin size at which
%       Amatrix(:,j) has been detected.
%       Unit_order  - unit order in Amatrix
%
%
%%  © 2020 Russo
%  for information please contact eleonora.russo@zi-mannheim.de
%%
if nargin<4 , display='raw'; end
if nargin == 4
    if ~strcmp(display,'raw') && ~strcmp(display,'ordunit') && ~strcmp(display,'clustered')
        fprintf('''raw'', ''ordunit'', ''clustered'' are the only acceptable style specifications\n')
        return
    end
end

nAss_final=size(As_across_bins,2);
Binvector=nan(1,nAss_final);

% Sparse-first construction to reduce peak RAM for large nneu x nAss matrices.
rows=[];
cols=[];
vals=[];
for i=1:nAss_final
    elems=double(As_across_bins{i}.elements(:));
    lags=double(As_across_bins{i}.lag(:));
    rows=[rows; elems]; %#ok<AGROW>
    cols=[cols; i*ones(numel(elems),1)]; %#ok<AGROW>
    vals=[vals; lags]; %#ok<AGROW>
    Binvector(i)=As_across_bins{i}.bin;
end
AAT_sparse=sparse(rows,cols,vals,nneu,nAss_final);
AAT=full(AAT_sparse);
AAT(AAT==0)=nan;

switch display
    case 'raw'
        Unit_order=1:nneu;
        As_order=1:length(As_across_bins);

    case 'ordunit'
        is_empty=all(isnan(AAT),2);
        idx_nan=find(is_empty);
        idx_activ=find(~is_empty);
        AAT=[AAT(~is_empty,:);AAT(is_empty,:)];
        Unit_order=[idx_activ;idx_nan];
        As_order=1:length(As_across_bins);

    case 'clustered'
        is_empty=all(isnan(AAT),2);
        idx_nan=find(is_empty);
        idx_activ=find(~is_empty);
        AAT_activ=AAT(~is_empty,:);

        A01=~isnan(AAT_activ);
        M_assemb=zeros(size(AAT_activ,1));
        for n=1:size(A01,2)
            aus=find(A01(:,n)==1);
            M_assemb(aus,aus)=M_assemb(aus,aus)+1;
        end

        M_assemb(M_assemb==0)=0.0001;
        d_assemb=ones(size(AAT_activ,1))./M_assemb;
        d_assemb=d_assemb-diag(diag(d_assemb));
        Q=linkage(d_assemb,'average');
        [~,~,perm]=dendrogram(Q,0);
        perm1=idx_activ(perm);

        D=pdist(double(~isnan(AAT_activ))');
        Z=squareform(D);
        Q2=linkage(Z,'average');
        [~,~,perm2]=dendrogram(Q2,0);

        AAT=AAT_activ(perm,perm2);
        AAT=[AAT; nan(length(idx_nan),size(AAT,2))];
        Unit_order=[perm1;idx_nan];
        Binvector=Binvector(perm2);
        As_order=perm2;
end

Amatrix=AAT;
binmat=log10(Binvector);
AAT(end+1,end+1)=0;
binmat(end+1,end+1)=0;

subplot(2,1,1)
h=pcolor(AAT);
set(h, 'EdgeColor', [0.8 0.8 0.8]);
caxis([0, max(AAT(:))])
colormap jet
ylabel('Unit #')
set(gca, 'XTick', []);
hcb=colorbar;
hcb.Label.String = 'Time lag \it l \rm (# bins)';
yy = 1:size(AAT,1);
set(gca,'YTick',yy(1)+0.5:yy(end),'Yticklabel',Unit_order)
set(gcf, 'Color', [1,1,1]);

subplot(2,1,2)
pcolor(binmat)
yy = 1:size(AAT,2);
set(gca,'XTick',yy(1)+0.5:2:yy(end)+1+0.05,'Xticklabel',As_order(1:2:end))
xlabel('Assembly #')
set(gca, 'YTick', []);
hC = colorbar;
if length(BinSizes)>1
    caxis([log10(min(BinSizes)), log10(max(BinSizes))]);
else
    caxis([log10(BinSizes-0.001), log10(BinSizes+0.001)]);
end
L=[0.001,0.002,0.003,0.004,0.005, 0.006,0.007,0.008,0.009,0.01,0.02,0.03,...
    0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,...
    1,2,3,4,5,6,7,8,9,10];
set(hC,'Ytick',log10(L),'YTicklabel',L);
set(hC,'Location','southoutside')
hC.Label.String = 'Temporal precision \Delta (sec)';