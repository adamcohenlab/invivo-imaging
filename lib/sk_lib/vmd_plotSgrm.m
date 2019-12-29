function [h] = vmd_plotSgrm(sgrm)
%plot spectrogram
%   to do: document alxis handles naming conventions
%       label plots
%       options for lin/log plotting of spectrogram
%       Simon kheifets 6/10/2018
Stot = sgrm.Stot;
fvecS=sgrm.fvec;
tvecS = sgrm.tvec;
fvecSLog=sgrm.fvecLog;
StotLog = sgrm.StotLog;
f = sgrm.fvm;
Svm = sgrm.Svm;

figure('position',[208         217        1437         685],'Color','w');
ah1 = subplot(2,3,[1 2 4 5]);
%imagesc(log10(Stot),'XData',tvecS,'YData',fvecS);
imagesc(tvecS,fvecS, 10*log10(Stot));
%imagesc(tvecS,fvecSLog, 10*log10(StotLog));
colorbar();
%set(gca,'CLim',[38 50])
title(sgrm.parent.label);
xlabel('time(s)');
ylabel('f(Hz)');

ah2 = subplot(2,3,3);
[h.xbin, h.ybin] = lowpassplot(f,Svm,300,'bintype','log','Axes',ah2);
xlim([.1 500]);
ylabel('PSD (counts^2/Hz)');
title('power density')

ah3 = subplot(2,3,6);
semilogx(f(f>10),cumsum(Svm(f>10)));
xlim([.1 500]);
title('cumulative power');
xlabel('frequency (Hz)');
ylabel('variance (counts^2)');

h.ah1 = ah1;
h.ah2 = ah2;
h.ah3 = ah3;

end

