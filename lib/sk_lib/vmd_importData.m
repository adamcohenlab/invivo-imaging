function [v] = vmd_importData(s)
%VMD_importData. Import data and create a vmd object based on settings 
% stored in struct s. create s using s = makeVMDsettings(...)
%   Currently written for files from upright in-vivo microscope
%   Simon kheifets 6/3/2008
    v = VMD();
    v.timestamp = datestr(now);
    [mov0, nrow, ncol]=readBinMov4(s.dir,s.file,s.transpose);
    v.mov=vm(mov0((1+s.crop(3)):(end-s.crop(4)),(1+s.crop(1)):(end-s.crop(2)),:));
    v.dt = s.dt;
    v.tvec = s.t0+s.dt*(1:v.mov.frames)-s.dt;
    v.meanimg = v.mov.mean;
    v.meantrace = v.mov.frameAverage;
    v.vartrace = var(double(v.mov(:,:)),0,1)'; %variance of each frame
    v.duration = v.tvec(end);
    if isempty(s.label)
        v.label = s.file;
    else
        v.label = s.label;
    end
    v.set = s;
    
    if s.calcvar
            v.varimg = var(single(v.mov.data),0,3);
    end
end

