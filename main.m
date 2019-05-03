labpath = '/n/cohen_lab/Lab';

addpath(fullfile(labpath,'Labmembers','Michael Xie','in vivo data processing','quivvia','lib'));
addpath(fullfile(labpath,'Computer Code','Image Processing'));
addpath(fullfile(labpath,'Computer Code','Image Processing','NoRMCorre-master'));
addpath(fullfile(labpath, 'Labmembers', 'Yoav Adam', 'Scripts', 'NoRMCorre-master'));

%%
% home = fullfile(labpath,'Labmembers/Linlin Fan/In vivo/IVOP14/2018-06-03_D24/Vol img/slice2/FOV21/185008_EW');
% home = fullfile(labpath,'Labmembers/Linlin Fan/In vivo/IVOP21/2018-12-02/slice1/FOV6/195745_PuffTE2');
% home = fullfile(labpath,'Labmembers/Linlin Fan/In vivo/IVOP21/2018-12-02/slice1/FOV14motion/204040_PuffTE2');
% home = fullfile(labpath,'Data and Analysis/Electrochromic Protein/somArchon optopatch in vivo/NMFcorrection_MEX/IVOP21/2019-01-08/slice2/FOV4/161916_PuffTE2');

% home = fullfile(labpath,'Labmembers','Linlin Fan','In vivo','IVOP14','2018-06-03_D24','Vol img','slice2','FOV21','184947_F');
home = '/n/cohen_lab/Lab/Labmembers/Yoav Adam/Data/In Vivo/PlaceCells/PC1R1/2018-11-16_PC1R1-S1/slice1/FOV3/131943_FreeRun_Dilas-8V_488-OD1.0-Mask0-Pos22';

output = fullfile(home,'PMD_output_MX');
if ~exist(output,'dir')
    mkdir(output)
end

%% extract motion traces into MAT file
reg_shifts = returnShifts(home);
save(fullfile(home,'reg_shifts.mat'),'reg_shifts');

%% denoising
mov_in = "movReg.bin";
detr_spacing = 5000;
row_blocks = 4;
col_blocks = 2;
stim_dir = [];...fullfile('/','matlab wvfm','PuffTE2','AOwaveforms.bin');
    
trunc_start = 1; % frame to start truncated movie
trunc_length = 6000; % length of truncated movie

%%
run_command = sprintf("cd denoise\n source setup.sh\n sbatch denoise.run ""%s"" ""%s"" ""%s"" %d %d %d %d %d ""%s""",...
    home, mov_in, output, detr_spacing, row_blocks, col_blocks,...
    trunc_start-1, trunc_length, stim_dir);

system(run_command);

%% motion correction
load(fullfile(home,'reg_shifts.mat'));
xShifts = reg_shifts(1,71:end);
yShifts = reg_shifts(2,71:end);
dX = xShifts - mean(xShifts);
dY = yShifts - mean(yShifts);
dXhp = dX - smooth(dX, 2000)';  % high pass filter
dYhp = dY - smooth(dY, 2000)';
dXs = smooth(dXhp, 5)';  % low-pass filter just to remove some jitter in the tracking.  Not sure if necessary
dYs = smooth(dYhp, 5)';


% mov = load(fullfile(output,'denoised.mat'));
mov = shiftdim(double(vm(fullfile(output,'denoised.tif'))),2);
[ySize, xSize, nFrames] = size(mov);
t = 1:nFrames;

% remove DMD refresh frames
intens = squeeze(mean(mean(mov)));
figure(880);
intensS = intens - smooth(intens, 100);
plot(intensS(100:end));title('pick threshold to remove frames');pause(1)
Rem=input('Pick thres to remove frames   ');
badFrames = find(intensS < Rem);
mov(:,:,badFrames) = mov(:,:,badFrames - 1);
intens1=squeeze(mean(mean(mov)));
figure(880);clf;hold on
plot(intens,'b')
plot(intens1,'r');

avgImg = mean(mov,3);
dmov = mov - avgImg;

dT = detr_spacing;
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


%% blood removal
figure(881); clf; moviefixsc(out4);
refimg = max(out4(:,:,1000:2000),[],3);

nframes = size(out4, 3);

figure(882); clf;
imshow(refimg, [], 'InitialMagnification', 'fit')
title('click to remove blood')
hold on;

inpoly = zeros(size(refimg));

[ysize, xsize] = size(refimg(:,:,1));
npts = 1;
colorindex = 0;
order = get(gca,'ColorOrder');
nroi = 1;
[x, y] = meshgrid(1:xsize, 1:ysize);
while(npts > 0)
    [xv, yv] = (getline(gca, 'closed'));
    if size(xv,1) < 3  % exit loop if only a line is drawn
        break
    end
    inpoly = inpoly + inpolygon(x,y,xv,yv);
    
    %draw the bounding polygons and label them
    currcolor = order(1+mod(colorindex,size(order,1)),:);
    plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
    text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
    
    colorindex = colorindex+1;
    roi_points{nroi} = [xv, yv];
    nroi = nroi + 1;
end

mov = out4.*repmat(inpoly==0, [1, 1, nframes]);
bloodmask = uint8(mean(mov,3) ~= 0);
saveastiff(bloodmask,fullfile(output,'bloodmask.tif'));

%% background selection
refimg = max(mov(:,:,1000:2000),[],3);
nframes = size(mov, 3);

figure(883); clf;
imshow(refimg, [], 'InitialMagnification', 'fit')
title('click to select background')
hold on;

inpoly = zeros(size(refimg));

[ysize, xsize] = size(refimg(:,:,1));
npts = 1;
colorindex = 0;
order = get(gca,'ColorOrder');
nroi = 1;
intens = [];
[x, y] = meshgrid(1:xsize, 1:ysize);
while(npts > 0)
    [xv, yv] = (getline(gca, 'closed'));
    if size(xv,1) < 3  % exit loop if only a line is drawn
        break
    end
    inpoly = inpoly + inpolygon(x,y,xv,yv);
    
    %draw the bounding polygons and label them
    currcolor = order(1+mod(colorindex,size(order,1)),:);
    plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
    text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
    
    colorindex = colorindex+1;
    roi_points{nroi} = [xv, yv];
    nroi = nroi + 1;
end

background = mov.*repmat(inpoly~=0, [1, 1, nframes]);
background = background - repmat(mean(background,3),[1 1 nframes]);
[U, S, V] = svds(double(reshape(background,[size(background,1)*size(background,2), nframes])),6);
ff = (V - mean(V,2));
fb = (U * S);
figure(884);stackplot(ff);
saveastiff(ff,fullfile(output,'ff.tif'))
saveastiff(fb,fullfile(output,'fb.tif'))
