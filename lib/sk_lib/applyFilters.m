function [D,V,extra] = applyFilters(Uvm,obj,tvec)
%apply filters (images stored in Uvm) to movie obj
%   Detailed explanation goes here
%   Simon kheifets 5/15/2018
totvar = sum(obj.data(:).^2);

Vnn = (Uvm(:,:)'*obj(:,:))'; %use the filters on the original movie
D = (sum(Vnn.^2,1))'; %amplitudes for PC components
V = Vnn./(ones(length(tvec),1)*sqrt(D')); %normalized time traces

%calculate autocorrelation for each time trace
tvecthis = tvec-tvec(1); tC = [-fliplr(tvecthis(2:end)) tvecthis];
Vcorr=zeros(length(tC),Uvm.frames);
for j = 1:size(V,2)
    v = V(:,j); 
    Vcorr(:,j) = xcorr(v-mean(v));
end
%stackplot(Vcorr(tC>0,:)); set(gca,'XScale','log');

%calculate PSD for each time trace
[Sv, f, df] = calc_PSD(V,hamming(length(tvec)),tvec(2)-tvec(1));

try
    Svbin = binfunc(Sv,'linlog',300,10);
    fbin = binfunc(f,'linlog',300,10);
catch
    Svbin = Sv;
    fbin = f;
end
% SvBinStack = Svbin.*(ones(length(fbin),1)*10.^-(1:nPCs));
% loglog(fbin,SvBinStack); xlim([10 500]);

extra.totvar = totvar;
extra.Sv = Sv;
extra.f = f;
extra.Svbin = Svbin;
extra.fbin = fbin;
extra.Vcorr = Vcorr;
extra.tC = tC;
end

