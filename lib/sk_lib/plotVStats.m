function [] = plotVStats(Dnew,Vnew,Exnew)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

nPCs = size(Vnew,2);

subplot(2,2,1)
stackplot(Vnew);
title('temporal components');

subplot(2,2,3)
Dc = [0; cumsum(Dnew)];
semilogy(0:1:length(Dnew),(Exnew.totvar-Dc)/Exnew.totvar,'o-');
title('fractional variance');

subplot(2,2,2)
semilogx(Exnew.tC, Exnew.Vcorr);
title('autocorrelation');

subplot(2,2,4)
SvBinStack = Exnew.Svbin.*(ones(length(Exnew.fbin),1)*10.^-(1:nPCs));
loglog(Exnew.fbin,SvBinStack); %xlim([10 500]);
title('Spectra');
end

