clc;clearvars;close all

addpath('/n/cohen_lab/Lab/Labmembers/Michael Xie/in vivo data processing/quivvia/lib');
addpath('/n/cohen_lab/Lab/Computer Code/Image Processing');
addpath('/n/cohen_lab/Lab/Computer Code/Image Processing/NoRMCorre-master')
labpath = labPath();
addpath(fullfile(labpath, 'Labmembers', 'Yoav Adam', 'Scripts', 'NoRMCorre-master'));

data_dir = 'data/';
%%
if exist(fullfile(data_dir,'reg_shifts.mat'),'file')
    regShifts = load(fullfile(data_dir,'reg_shifts.mat'));
    xShifts = regShifts.reg_shifts(1,:);
    yShifts = regShifts.reg_shifts(2,:);
    nFrames = size(xShifts,2);
    dX = xShifts - mean(xShifts);
    dY = yShifts - mean(yShifts);
    dXhp = dX - smooth(dX, 2000)';  % high pass filter
    dYhp = dY - smooth(dY, 2000)';
    dXs = smooth(dXhp, 5)';  % low-pass filter just to remove some jitter in the tracking.  Not sure if necessary
    dYs = smooth(dYhp, 5)';
    
    
    mov = vm(fullfile(data_dir,'denoised.tif'));
    mov = shiftdim(double(mov), 2);
    [ySize, xSize, nFrames] = size(mov);
    t = 1:nFrames;
    
    avgImg = mean(mov,3);
    dmov = mov - avgImg;
    
    dT = 5000;
    % First column is the start of each epoch, second column is the end
    bdry = [(1:dT:nFrames)', [(dT:dT:nFrames) nFrames]'];
    nepoch = size(bdry, 1);
    out4 = zeros(size(mov));
    for j = 1:nepoch;
        tau = bdry(j,1):bdry(j,2);
        [out4(:,:,tau), ~] = SeeResiduals(dmov(:,:,tau), [dXs(tau); dYs(tau); dXs(tau).^2; dYs(tau).^2; dXs(tau) .* dYs(tau)], 1);
    end;
   
    saveastiff(out4,fullfile(data_dir,'motion_corrected.tif'));
end

exit;