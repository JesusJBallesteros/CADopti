function [As_across_bins,As_across_bins_index,assembly,Assemblies_all_orders]=CADopti(spM,MaxLags,BinSizes,ref_lag,alph,No_th,O_th,bytelimit)
% this function returns cell assemblies detected in spM spike matrix binned 
% at a temporal resolution specified in 'BinSizes' vector and testing for all 
% lags between '-MaxLags(i)' and 'MaxLags(i)'
%
% USAGE: [assembly]=Main_assemblies_detection(spM, MaxLags, BinSizes, ref_lag, alph, Dc, No_th, O_th, bytelimit)
%
% ARGUMENTS:
% spM     := matrix with population spike trains; each row is the spike train (time stamps, not binned) relative to a unit. 
% BinSizes:= vector of bin sizes to be tested;
% MaxLags:= vector of maximal lags to be tested. For a binning dimension of BinSizes(i) the program will test all pairs configurations with a time shift between -MaxLags(i) and MaxLags(i);
% (optional) ref_lag      := reference lag. Default value 2
% (optional) alph      := alpha level. Default value 0.05
% (optional) No_th     := minimal number of occurrences required for an assembly (all assemblies, even if significant, with fewer occurrences than No_th are discarded). Default value 0.
% (optional) O_th      := maximal assembly order (the algorithm will return assemblies of composed by maximum O_th elements).
% (optional) bytelimit := maximal size (in bytes) allocated for all assembly structures detected with a bin dimension. When the size limit is reached the algorithm stops adding new units. 
%
% RETURNS:
% assembly - structure containing assembly information:
%     assembly.parameters       - parameters used to run Main_assemblies_detection
%     assembly.bin{i} contains information about assemblies detected with 
%                     'BinSizes(i)' bin size tested for all lags between 
%                     '-MaxLags(i)' and 'MaxLags(i)'
% 
%        assembly.bin{i}.bin_edges - bin edges (common to all assemblies in assembly.bin{i})
%        assembly.bin{i}.n{j} information about the j-th assembly detected with BinSizes(i) bin size 
%                   elements: vector of units taking part to the assembly (unit order correspond to the agglomeration order)
%                        lag: vector of time lags. '.lag(z)' is the activation delay between .elements(1) and .elements(z+1)
%                         pr: vector of pvalues. '.pr(z)' is the pvalue of the statistical test between performed adding .elements(z+1) to the structure .elements(1:z)
%                       Time: assembly activation time. It reports how many times the complete assembly activates in that bin. .Time always refers to the activation of the first listed assembly element (.elements(1)), that doesn't necessarily corresponds to the first unit firing.
%               Noccurrences: number of assembly occurrence. '.Noccurrences(z)' is the occurrence number of the structure composed by the units .elements(1:z+1) 
%
% As_across_bins         - structure containing assembly information (exactly same information contained in "assembly" but collected across different temporal resolutions)
% As_across_bins_index   - information to link assemblies in "As_across_bins"
%                          back to the structure "assembly":
%                          assembly As_across_bins{i} is contained in assembly.bin{As_across_bins_index{i}(1)}.n{As_across_bins_index{i}(2)}.
%
%  � 2020 Russo
%  for information please contact eleonora.russo@zi-mannheim.de

if nargin<4 || isempty(ref_lag), ref_lag=2; end
if nargin<5 || isempty(alph), alph=0.05; end
if nargin<6 || isempty(No_th), No_th=0; end
if nargin<7 || isempty(O_th), O_th=Inf; end
if nargin<8 || isempty(bytelimit), bytelimit=Inf; end

nneu=size(spM,1);
nBins=length(BinSizes);
testit=ones(1,nBins);
binM=cell(1,nBins);
lag_tests=(2*MaxLags+1);
number_tests=0;
bin_bytes_used=zeros(1,nBins);

% matrix binning at all bins
for gg=1:nBins
    int=BinSizes(gg);
    tb=min(spM(:)):int:max(spM(:));

    binned_counts=zeros(nneu,length(tb)-1,'uint16');
    number_tests=number_tests+nneu*(nneu-1)*lag_tests(gg)/2;
    for n=1:nneu
        [binned_counts(n,:),~]=histcounts(spM(n,:),tb);
    end

    % RAM-aware storage: logical when all counts are binary, uint8 otherwise.
    if all(binned_counts(:)<=1)
        binM{gg}=logical(binned_counts);
    else
        binM{gg}=uint8(min(binned_counts,uint16(intmax('uint8'))));
    end

    assembly.bin{gg}.n=[];
    assembly.bin{gg}.bin_edges=tb;
    if size(binM{gg},2)-MaxLags(gg)<100
        fprintf('Warning: testing bin size=%f. The time series is too short, consider taking a longer portion of spike train or diminish the bin size to be tested \n', int);
        testit(gg)=0;
    end
end

fprintf('order 1\n')
clear Assemblies_all_orders
O=1;
Dc=100;

npairs=nneu*(nneu-1)/2;
pair_w1=zeros(npairs,1,'uint32');
pair_w2=zeros(npairs,1,'uint32');
k=1;
for w1=1:nneu
    for w2=w1+1:nneu
        pair_w1(k)=w1;
        pair_w2(k)=w2;
        k=k+1;
    end
end

p_values_pairs=nan(npairs,nBins);
assembly_selected_xy=cell(1,npairs);

parfor pair_idx=1:npairs
    w1=pair_w1(pair_idx);
    w2=pair_w2(pair_idx);
    assemblybin=cell(1,nBins);
    p_by_bin=nan(1,nBins);
    for gg=1:nBins
        assemblybin{gg}=FindAssemblies_recursive_prepruned([double(binM{gg}(w1,:));double(binM{gg}(w2,:))],w1,w2,MaxLags(gg),Dc,ref_lag);
        p_values_pairs(pair_idx,gg)=assemblybin{gg}.pr(end);
        assemblybin{gg}.bin=BinSizes(gg);
        assemblybin{gg}.bin_idx=gg;
        p_by_bin(gg)=assemblybin{gg}.pr;
    end
    [~, b]=min(p_by_bin);
    assembly_selected_xy{pair_idx}=compact_assembly_struct(assemblybin{b});
end

p_values=reshape(p_values_pairs.',[],1);
assembly_selected=cell(1,npairs);
keep_idx=1;
for pair_idx=1:npairs
    candidate=assembly_selected_xy{pair_idx};
    if candidate.Noccurrences(end) < No_th
        continue;
    end
    candidate_bytes=estimate_struct_bytes(candidate);
    bin_idx=double(candidate.bin_idx);
    if bin_bytes_used(bin_idx)+candidate_bytes > bytelimit
        continue;
    end
    bin_bytes_used(bin_idx)=bin_bytes_used(bin_idx)+candidate_bytes;
    assembly_selected{keep_idx}=candidate;
    keep_idx=keep_idx+1;
end
assembly_selected(keep_idx:end)=[];

if ~isempty(assembly_selected)

x=1:length(p_values);
p_values=sort(p_values);
p_values_alpha=alph./(number_tests+1-x);
aus=find((p_values'-p_values_alpha)<0);
if isempty(aus)
    HBcorrected_p=0;
else
    HBcorrected_p=p_values(aus(end));
end
ANfo=zeros(nneu,nneu,'logical');

keep_mask=true(1,length(assembly_selected));
for oo=1:length(assembly_selected)
    if assembly_selected{oo}.pr(end)>HBcorrected_p
        bin_idx=double(assembly_selected{oo}.bin_idx);
        bin_bytes_used(bin_idx)=max(0,bin_bytes_used(bin_idx)-estimate_struct_bytes(assembly_selected{oo}));
        keep_mask(oo)=false;
    else
        ANfo(double(assembly_selected{oo}.elements(1)),double(assembly_selected{oo}.elements(2)))=true;
    end
end
assembly_selected=assembly_selected(keep_mask);
Assemblies_all_orders{O}=assembly_selected;

Oincrement=1;
flush_chunk=512;
while Oincrement && O<(O_th-1)
    O=O+1;
    fprintf('order %d\n',O)
    Oincrement=0;
    assembly_selected_aus=cell(1,flush_chunk);
    xx=1;
    p_values_new=cell(1,size(assembly_selected,2));

    for w1=1:size(assembly_selected,2)
        ggg=double(assembly_selected{w1}.bin_idx);

        w1_elements=double(assembly_selected{w1}.elements);
        [~, w2_to_test]=find(ANfo(w1_elements,:)==1);
        w2_to_test(ismember(w2_to_test,w1_elements))=[];
        w2_to_test=unique(w2_to_test);
        local_pvals=nan(1,length(w2_to_test));

        for ww2=1:length(w2_to_test)
            w2=w2_to_test(ww2);
            spikeTrain2=double(binM{ggg}(w2,:))';
            assemblybin_aus=TestPair_ref(assembly_selected{w1},spikeTrain2,w2,MaxLags(ggg),Dc,ref_lag);
            local_pvals(ww2)=assemblybin_aus.pr(end);
            number_tests=number_tests+lag_tests(ggg);
            if assemblybin_aus.pr(end)<HBcorrected_p && assemblybin_aus.Noccurrences(end)>=No_th
                assemblybin_aus.bin=BinSizes(ggg);
                assemblybin_aus.bin_idx=ggg;
                assemblybin_aus=compact_assembly_struct(assemblybin_aus);
                candidate_bytes=estimate_struct_bytes(assemblybin_aus);
                if bin_bytes_used(ggg)+candidate_bytes<=bytelimit
                    if xx>numel(assembly_selected_aus)
                        assembly_selected_aus{xx+flush_chunk-1}=[];
                    end
                    assembly_selected_aus{xx}=assemblybin_aus;
                    bin_bytes_used(ggg)=bin_bytes_used(ggg)+candidate_bytes;
                    xx=xx+1;
                    Oincrement=1;
                end
            end
        end
        p_values_new{w1}=local_pvals;
    end

    if xx>1
        assembly_selected_aus(xx:end)=[];
    else
        assembly_selected_aus={};
    end

    if ~isempty(p_values_new)
        for pv_i=1:length(p_values_new)
            if ~isempty(p_values_new{pv_i})
                p_values=[p_values; p_values_new{pv_i}(:)]; %#ok<AGROW>
            end
        end
    end

    if Oincrement
        [assembly_selected,bin_bytes_used]=prune_same_size_assemblies(assembly_selected_aus,bin_bytes_used);
        Assemblies_all_orders{O}=assembly_selected;
    end

    x=1:length(p_values);
    p_values=sort(p_values);
    p_values_alpha=alph./(number_tests+1-x);
    aus=find((p_values'-p_values_alpha)<0);
    if isempty(aus)
        HBcorrected_p=0;
    else
        HBcorrected_p=p_values(aus(end));
    end
end

x=1:length(p_values);
p_values=sort(p_values)';
p_values_alpha=alph./(number_tests+1-x);
aus=find((p_values-p_values_alpha)<0);
if isempty(aus)
    HBcorrected_p=0;
else
    HBcorrected_p=p_values(aus(end));
end

for o=1:length(Assemblies_all_orders)
    keep_mask=true(1,length(Assemblies_all_orders{o}));
    for oo=1:length(Assemblies_all_orders{o})
        if Assemblies_all_orders{o}{oo}.pr(end)>HBcorrected_p
            bin_idx=double(Assemblies_all_orders{o}{oo}.bin_idx);
            bin_bytes_used(bin_idx)=max(0,bin_bytes_used(bin_idx)-estimate_struct_bytes(Assemblies_all_orders{o}{oo}));
            keep_mask(oo)=false;
        end
    end
    Assemblies_all_orders{o}=Assemblies_all_orders{o}(keep_mask);
end

o=length(Assemblies_all_orders);
x=1;
Element_template=cell(1,max(1,length(Assemblies_all_orders{o})));
for oo=1:length(Assemblies_all_orders{o})
    Element_template{x}=Assemblies_all_orders{o}{oo}.elements;
    x=x+1;
end

for o=length(Assemblies_all_orders)-1:-1:1
    keep_mask=true(1,length(Assemblies_all_orders{o}));
    for oo=1:length(Assemblies_all_orders{o})
        found=0;
        for ooo=1:(x-1)
            if ismember(double(Assemblies_all_orders{o}{oo}.elements),double(Element_template{ooo}))
                bin_idx=double(Assemblies_all_orders{o}{oo}.bin_idx);
                bin_bytes_used(bin_idx)=max(0,bin_bytes_used(bin_idx)-estimate_struct_bytes(Assemblies_all_orders{o}{oo}));
                keep_mask(oo)=false;
                found=1;
                break
            end
        end
        if found==0
            if x>numel(Element_template)
                Element_template{x+256}=[];
            end
            Element_template{x}=Assemblies_all_orders{o}{oo}.elements;
            x=x+1;
        end
    end
    Assemblies_all_orders{o}=Assemblies_all_orders{o}(keep_mask);
end

for o=length(Assemblies_all_orders):-1:1
    for oo=length(Assemblies_all_orders{o}):-1:1
        bx=double(Assemblies_all_orders{o}{oo}.bin_idx);
        assembly.bin{bx}.n=[assembly.bin{bx}.n,Assemblies_all_orders{o}(oo)]; %#ok<AGROW>
    end
end
for gg=length(BinSizes):-1:1
    if isempty(assembly.bin{gg}.n)
        assembly.bin{gg}=[];
    end
end

fprintf('\n');
else
    for gg=length(BinSizes):-1:1
        assembly.bin{gg}=[];
        Assemblies_all_orders=[];
    end
end

assembly.parameters.alph=alph;
assembly.parameters.Dc=Dc;
assembly.parameters.No_th=No_th;
assembly.parameters.O_th=O_th;
assembly.parameters.bytelimit=bytelimit;
assembly.parameters.ref_lag=ref_lag;

[As_across_bins,As_across_bins_index]=assemblies_across_bins(assembly,BinSizes);

end

function [assembly_out,bin_bytes_used]=prune_same_size_assemblies(assembly_in,bin_bytes_used)
if isempty(assembly_in)
    assembly_out={};
    return
end
na=length(assembly_in);
nelement=size(assembly_in{1}.elements,2);
selection=nan(na,nelement+1);
assembly_out=cell(1,na);
nns=1;
for i=1:na
    elem=sort(double(assembly_in{i}.elements));
    [ism,indx]=ismember(elem,selection(:,1:nelement),'rows');
    if ~ism
        assembly_out{nns}=assembly_in{i};
        selection(nns,1:nelement)=elem;
        selection(nns,nelement+1)=assembly_in{i}.pr(end);
        nns=nns+1;
    else
        if selection(indx,nelement+1)>assembly_in{i}.pr(end)
            bin_idx_old=double(assembly_out{indx}.bin_idx);
            bin_bytes_used(bin_idx_old)=max(0,bin_bytes_used(bin_idx_old)-estimate_struct_bytes(assembly_out{indx}));
            assembly_out{indx}=assembly_in{i};
            selection(indx,nelement+1)=assembly_in{i}.pr(end);
        else
            bin_idx_drop=double(assembly_in{i}.bin_idx);
            bin_bytes_used(bin_idx_drop)=max(0,bin_bytes_used(bin_idx_drop)-estimate_struct_bytes(assembly_in{i}));
        end
    end
end
assembly_out(nns:end)=[];
end

function assem=compact_assembly_struct(assem)
if isfield(assem,'elements')
    assem.elements=uint32(assem.elements);
end
if isfield(assem,'lag')
    assem.lag=int16(assem.lag);
end
if isfield(assem,'Noccurrences')
    assem.Noccurrences=uint32(max(0,assem.Noccurrences));
end
if isfield(assem,'bin_idx')
    assem.bin_idx=uint16(assem.bin_idx);
end
if isfield(assem,'Time') && ~isa(assem.Time,'uint16')
    assem.Time=uint16(max(0,assem.Time));
end
end

function out_bytes=estimate_struct_bytes(in_struct)
% Fast approximate byte estimate to avoid frequent costly WHOS calls.
out_bytes=0;
if isfield(in_struct,'elements')
    out_bytes=out_bytes+4*numel(in_struct.elements);
end
if isfield(in_struct,'lag')
    out_bytes=out_bytes+2*numel(in_struct.lag);
end
if isfield(in_struct,'pr')
    out_bytes=out_bytes+8*numel(in_struct.pr);
end
if isfield(in_struct,'Time')
    out_bytes=out_bytes+2*numel(in_struct.Time);
end
if isfield(in_struct,'Noccurrences')
    out_bytes=out_bytes+4*numel(in_struct.Noccurrences);
end
out_bytes=out_bytes+128;
end
