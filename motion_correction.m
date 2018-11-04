clc;clearvars;close all

% addpath('/n/cohen_lab/Lab/Labmembers/Michael Xie/in vivo data processing/quivvia/lib');
% addpath('/n/cohen_lab/Lab/Computer Code/Image Processing');
% addpath('/n/cohen_lab/Lab/Computer Code/Image Processing/NoRMCorre-master')
% labpath = labPath();
% addpath(fullfile(labpath, 'Labmembers', 'Yoav Adam', 'Scripts', 'NoRMCorre-master'));

uprFOVs = [...
    29,2017,09,14,5,1,4;...#1
    29,2017,09,14,5,1,6;...#2
    32,2017,08,17,2,1,1;...#3
    32,2017,08,17,2,1,2;...#4
    32,2017,08,17,2,1,3;...#5
    32,2017,08,17,2,1,4;...#6
    32,2017,08,17,2,1,6;...#7
    32,2017,09,14,4,1,2;...#8
    32,2017,09,14,4,1,4;...#9
    38,2017,09,20,1,1,1;...#10
    38,2017,09,20,1,1,2;...#11
    38,2017,09,20,1,1,3;...#12
    38,2017,09,20,1,1,4;...#13
    38,2017,09,20,1,1,5;...#14
    38,2017,10,04,2,1,1;...#15
    38,2017,10,04,2,1,2;...#16
    38,2017,10,04,2,1,3;...#17
    38,2017,10,04,2,1,4;...#18
    38,2017,10,04,2,1,5;...#19
    38,2017,10,04,2,1,6;...#20
    38,2017,10,04,2,1,7;...#21
    39,2017,09,21,1,1,1;...#22
    39,2017,09,21,1,1,2;...#23
    39,2017,09,21,1,1,3;...#24
    39,2017,09,21,1,1,4;...#25
    39,2017,09,21,1,1,5;...#26
    39,2017,09,21,1,1,6;...#27
    39,2017,09,21,1,1,7;...#28
    39,2017,10,04,2,1,1;...#29
    39,2017,10,04,2,1,2;...#30
    39,2017,10,04,2,1,3;...#31
    39,2017,10,04,2,1,4;...#32
    39,2017,10,04,2,1,5;...#33
    39,2017,10,04,2,1,6;...#34
    32,2017,08,17,2,1,1;...#35
    48,2018,06,21,7,1,1;...#36
    48,2018,06,21,7,1,2;...#37
    48,2018,06,21,7,1,3;...#38
    48,2018,06,21,7,1,4;...#39
    48,2018,06,21,7,1,5;...#40
    48,2018,06,21,7,1,7;...#41
    48,2018,06,21,7,1,8]; %#42 %ivq# ,y,m,d, session,slice,fov,

fovpaths = fovPath(uprFOVs);
regalbasepath = './';
data_dirs = cellfun(@(x) fullfile(regalbasepath,x),fovpaths,'UniformOutput',false);
nMovs = size(data_dirs,2);
%%
% parpool('local', 2);
parfor_progress(nMovs);
for idx = 1:nMovs
    if exist(fullfile(data_dirs{idx},'reg_shifts.mat'),'file')
        regShifts = load(fullfile(data_dirs{idx},'reg_shifts.mat'));
        xShifts = regShifts.reg_shifts(1,:);
        yShifts = regShifts.reg_shifts(2,:);
        nFrames = size(xShifts,2);
        dX = xShifts - mean(xShifts);
        dY = yShifts - mean(yShifts);
        dXhp = dX - smooth(dX, 2000)';  % high pass filter
        dYhp = dY - smooth(dY, 2000)';
        dXs = smooth(dXhp, 5)';  % low-pass filter just to remove some jitter in the tracking.  Not sure if necessary
        dYs = smooth(dYhp, 5)';


        mov = vm(fullfile(data_dirs{idx},'denoised.tif'));
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
            [out4(:,:,tau), ~] = SeeResiduals(dmov(:,:,tau), [dXs(tau); dYs(tau); dXs(tau).^2; dYs(tau).^2], 1);
        end;

        saveastiff(out4,fullfile(data_dirs{idx},'motion_corrected.tif'));
    %     dxpsd(:,idx) = fft(dX).*conj(fft(dX));
    %     dypsd(:,idx) = fft(dY).*conj(fft(dY));
    %     freq(:,idx) = (0:(nFrames-1))*1000/nFrames;
    end
    parfor_progress;
end
parfor_progress(0);

% figure(5);
% subplot(121);
% fill([freq(2:end);flipud(freq(2:end))],[max(dxpsd_m(2:end)-0.2*dxpsd_sd(2:end),0.1);flipud(dxpsd_m(2:end)+0.2*dxpsd_sd(2:end))],[.9 .9 .9],'linestyle','none');
% line(freq,dxpsd_m);
% set(gca, 'XScale', 'log', 'YScale','log')
% xlim([1 500])
% xlabel('Frequency (Hz)')
% ylabel('Power (A.U.)')
% title('Avg. power spectrum of x-position fluctuations')
% subplot(122);
% fill([freq(2:end);flipud(freq(2:end))],[max(dypsd_m(2:end)-0.2*dypsd_sd(2:end),0.1);flipud(dypsd_m(2:end)+0.2*dypsd_sd(2:end))],[.9 .9 .9],'linestyle','none');
% line(freq,dypsd_m);
% set(gca, 'XScale', 'log', 'YScale','log')
% xlim([1 500])
% xlabel('Frequency (Hz)')
% ylabel('Power (A.U.)')
% title('Avg. power spectrum of y-position fluctuations')
