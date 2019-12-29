function [icsFilter,icsFootprint,icsTime,nICs,mixmat,sepmat,extra] = vmICA(Uvm,D,V,nICs,alph,tvec)
%vmICA do ICA on a factorized movie M=Uvm*D*V' giving (at most) nICs
%components. Does the spatiotemporal ICA, with 0<alph<1 determining the
%ratio of spatial/temporal contribution to independence
%   Detailed explanation goes here

comSig = [(1-alph)*Uvm(:,:)*diag(sqrt(D)); alph*V*diag(sqrt(D))];
[icsST, mixmat, sepmat] = sorted_ica2(comSig',nICs);
icsST = icsST';
nICs = size(icsST,2);
%icUvm = vm(icsST(1:size(Uvm(:,:),1),:),[Uvm.rows Uvm.cols]);
%stackImshow(icUvm);
vT = V';
%vnewT = Vnew';
icsTime = sepmat*vT;
%icsTimeNew = sepmat*vnewT(p_keep,:);

icsFilter = vm((Uvm(:,:))*sepmat',[Uvm.rows,Uvm.cols]); %filter
icsFootprint = vm((Uvm(:,:))*mixmat,[Uvm.rows,Uvm.cols]); %footprint

%statistical properties of icsTime components:

V=icsTime;
%calculate autocorrelation for each time trace
tvecthis = tvec-tvec(1); tC = [-fliplr(tvecthis(2:end)) tvecthis];
Vcorr=zeros(length(tC),size(V,1));
for j = 1:size(V,1)
    v = V(j,:); 
    Vcorr(:,j) = xcorr(v-mean(v));
end
%stackplot(Vcorr(tC>0,:)); set(gca,'XScale','log');

%calculate PSD for each time trace
[Sv, f, df] = calc_PSD(V',hamming(length(tvec)),tvec(2)-tvec(1));
try
    Svbin = binfunc(Sv,'linlog',300,10);
    fbin = binfunc(f,'linlog',300,10);
catch
    Svbin = Sv;
    fbin = f;
end
% SvBinStack = Svbin.*(ones(length(fbin),1)*10.^-(1:nPCs));
% loglog(fbin,SvBinStack); xlim([10 500]);

extra.totvar = 0;
extra.Sv = Sv;
extra.f = f;
extra.Svbin = Svbin;
extra.fbin = fbin;
extra.Vcorr = Vcorr;
extra.tC = tC;

end

