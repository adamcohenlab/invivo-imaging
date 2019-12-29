function [objout] = rmRowNoise(obj,tvec)
%RMROWNOISE removes noise that is coherent across each row
%   obj is a vm movie
%   tvec is a vector map of time points
%   Simon Kheifets 6/10/2018

dt = tvec(2)-tvec(1);
nlp = 40;
rowavg = squeeze(mean(obj.data,2));
nPCs= 20;
whichPC = 2:14; %which PCs to keep for recreating noise...

ravm = vm(rowavg,[size(rowavg,1),1]);
ravm = ravm-ravm.mean;
ravmLP = vm((filter(ones(1,nlp)/nlp,1,(ravm(:,:).data)'))',[ravm.rows,ravm.cols]);
ravmHP = ravm-ravmLP;

[Uvm,D,V,Ex] = vmPCA(ravmHP,nPCs,tvec); %do PCA
figure;
plotVStats(D,V,Ex)


noisecol = Uvm(:,whichPC)*sqrt(diag(D(whichPC)))*(V(:,whichPC))';
noisecolvm = vm(noisecol,[size(rowavg,1),1]);

noisecol = reshape(noisecol,obj.rows,1,obj.frames);
noisemov = vm(repmat(noisecol,1,obj.cols,1));
objout = obj-noisemov;
end

