function [Uvm,D,V,extra] = vmPCA(obj,nPCs,tvec,varargin)
%VMPCA  PCA filtering
%   Detailed explanation goes here
%szM = size(obj(:,:));
if nargin == 4
    doplot = varargin{1};
else
    doplot = 1;
end
[U,S,V] = svds(obj(:,:),nPCs);
Uvm = vm(U,[obj.rows,obj.cols]); %store PC images as a vm
totvar = sum(obj.data(:).^2);
D = (diag(S)).^2; 

%calculate autocorrelation for each time trace
tvecthis = tvec-tvec(1); tC = [-fliplr(tvecthis(2:end)) tvecthis];
Vcorr=zeros(length(tC),nPCs);
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

%display PC images, time traces, correlation functions, and PSD
if doplot
    figure;
    try
    [handles] = plotPCAresults(Uvm,V,D,totvar,'PCA');
    end
end

extra.totvar = totvar;
extra.Sv = Sv;
extra.f = f;
extra.Svbin = Svbin;
extra.fbin = fbin;
extra.Vcorr = Vcorr;
extra.tC = tC;
%% todo
%include spectra and correlograms in plotPCA
%manually select modes?
%
end

