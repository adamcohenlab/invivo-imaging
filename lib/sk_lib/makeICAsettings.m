function [s] = makeICAsettings(varargin)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
%   Simon kheifets 6/3/2018
p = inputParser;
p.KeepUnmatched = 1; %don't give error when unused paramters appera

addParameter(p,'NIC',[]);
addParameter(p,'PKeep',[]);
addParameter(p,'Alpha',[]);


parse(p,varargin{:});
v2struct(p.Results);

s.nic = NIC;
s.pkeep = PKeep;
s.alpha = Alpha;
end

