function [sgrm] = vmd_calcSgrm(obj)
%VMD_calcSgrm... obj is a vmd object
%   Simon kheifets 6/10/2018

rowavg = squeeze(mean(obj.mov.data,2));
[Sv, f, df] = calc_PSD(rowavg',hamming(length(rowavg)),obj.dt);
Svm = mean(Sv,2);
% loglog(f+df,Svm);
% semilogx(f(f>10),cumsum(Svm(f>10)));
% xlim([10 500]);
%spectrogram of rowavg

nsc = floor(length(rowavg)/10);
nov = floor(nsc/2);
nff = max(256,2^nextpow2(nsc));


for i = 1:size(rowavg,1)
    x = rowavg(i,:);
    S = spectrogram(x,hamming(nff),nov,nff,'MinThreshold',0);
    if i ==1
        Stot = abs(S);
    else
        Stot = Stot+abs(S);
    end
end
fvecS = linspace(0,500,size(S,1)+1);
fvecS = fvecS(1:end-1);
tvecS = linspace(obj.tvec(1),obj.tvec(end),size(S,2));

%log binning for log plotting
for i = 1:size(Stot,2)
    StotLog(:,i)=binfunc((Stot(:,i)),'log',100,1);
end
fvecSLog = binfunc(fvecS(:,i),'log',100,1);

sgrm.Stot = Stot;
sgrm.fvec=fvecS;
sgrm.tvec = tvecS;
sgrm.fvecLog=fvecSLog;
sgrm.StotLog = StotLog;
sgrm.parent = obj;
sgrm.Svm = Svm;
sgrm.fvm = f;
obj.sgrm = sgrm; %attatch spectrogram information to vmd object
end

