function [s] = makeVMDsettings(varargin)
%Is it possible to set some private settings of s here which can thus be
%defined in the code but not visible when s is listed.
%example, list of supported platforms
%version number
%   Detailed explanation goes here
%   Simon kheifets 6/3/2008
p = inputParser;
p.KeepUnmatched = 1; %don't give error when unused paramters appera
addParameter(p,'Dir',[]);
addParameter(p,'File',[]);
addParameter(p,'Transpose',1);
addParameter(p,'DT',1e-3);
addParameter(p,'T0',0);
addParameter(p,'CalcVar',1); %calculate variance during loading?
addParameter(p,'Label',[]);
%addParameter(p,'Platform','Windows'); %{'Windows'; 'Unix'}
addParameter(p,'Crop',[0 0 0 0]); %[leftpx rightpx bottompx toppx]
addParameter(p,'PixelSize',1);          %size of one pixel in pixelunits
addParameter(p,'PixelUnits', 'camera px');% 
parse(p,varargin{:});
v2struct(p.Results);

s.dir = Dir;
s.file = File;
s.transpose = Transpose;
s.dt = DT;
s.t0 = T0;
s.calcvar = CalcVar;
s.label = Label;
s.crop = Crop;
s.pixelsize = PixelSize;
s.pixelunits = PixelUnits;
    
end

