function [s] = makeFVMsettings(varargin)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
%   Simon kheifets 6/3/2018
p = inputParser;
p.KeepUnmatched = 1; %don't give error when unused paramters appera

addParameter(p,'TLim',[]);
addParameter(p,'Label','');
addParameter(p,'CropMotion',1); %crop motion edges
addParameter(p,'DoPBLC',1); %do photobleach correction?
addParameter(p,'RemoveRowNoise',1); %remove row noise?
addParameter(p,'NBin',3); %binning number
addParameter(p,'FixDMDFrames',1);
addParameter(p,'HPSamples',20); %time constant (in samples) of highpass filter
addParameter(p,'RemoveBloodflow',1); %try to mask pixels w bloodflow
addParameter(p,'TBlood',2.5); %threshold for variance of pixels to remove


parse(p,varargin{:});
v2struct(p.Results);

s.tlim = TLim;
s.label = Label;
s.cropmotion = CropMotion;
s.dopblc = DoPBLC;
s.removerownoise = RemoveRowNoise;
s.nbin = NBin;
s.fixdmdframes = FixDMDFrames;
s.hpsamples = HPSamples;
s.rmbloodflow = RemoveBloodflow;
s.tblood = TBlood;
    
end
