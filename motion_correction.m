% clc;clearvars;close all

labpath = '/n/cohen_lab/Lab';

addpath(fullfile(labpath,'Labmembers','Michael Xie','in vivo data processing','quivvia','lib'));
addpath(fullfile(labpath,'Computer Code','Image Processing'));
addpath(fullfile(labpath,'Computer Code','Image Processing','NoRMCorre-master'));
addpath(fullfile(labpath, 'Labmembers', 'Yoav Adam', 'Scripts', 'NoRMCorre-master'));

%%
if exist(fullfile(home,'reg_shifts.mat'),'file')
    load(fullfile(home,'reg_shifts.mat'));
    xShifts = reg_shifts(1,71:end);
    yShifts = reg_shifts(2,71:end);
    dX = xShifts - mean(xShifts);
    dY = yShifts - mean(yShifts);
    dXhp = dX - smooth(dX, 2000)';  % high pass filter
    dYhp = dY - smooth(dY, 2000)';
    dXs = smooth(dXhp, 5)';  % low-pass filter just to remove some jitter in the tracking.  Not sure if necessary
    dYs = smooth(dYhp, 5)';
    
    tic;
    while(~exist(fullfile(output,'PMD_residual.tif'),'file'))
        pause(30);
        if(toc > 24 * 60)
            exit;
        end
    end
    
    tic;
    mov = shiftdim(double(vm(fullfile(output,'denoised.tif'))),2);
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
    
    motion_corrected = fullfile(output,'motion_corrected.bin');
    fid = fopen(motion_corrected,'w');
    fwrite(fid,single(out4),'float32');
    fclose(fid);
    toc;
end

exit;