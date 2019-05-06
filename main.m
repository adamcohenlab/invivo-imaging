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

%% denoising parameters
mov_in = "movReg.bin";
detr_spacing = 5000;
row_blocks = 4;
col_blocks = 2;
stim_dir = [];...fullfile('/','matlab wvfm','PuffTE2','AOwaveforms.bin');
    
trunc_start = 1; % frame to start truncated movie
trunc_length = 6000; % length of truncated movie

%% denoising and motion correction

run_command = sprintf("cd denoise\n source setup.sh\n sbatch denoise.run ""%s"" ""%s"" ""%s"" %d %d %d %d %d ""%s""",...
    home, mov_in, output, detr_spacing, row_blocks, col_blocks,...
    trunc_start-1, trunc_length, stim_dir);

system(run_command);

%% blood removal

noise_im = loadtiff(fullfile(output,'Sn_image.tif'));
[ysize, xsize] = size(noise_im);

out4 = single(vm(fullfile(output,'motion_corrected.bin'),ysize,xsize));

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
