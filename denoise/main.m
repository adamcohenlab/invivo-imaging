%%
addpath(genpath(fullfile('..','lib'));

%%
home = fullfile('..','demo_data');

output = fullfile(home,'output');

if ~exist(output,'dir')
    mkdir(output)
end

%% NoRMCorre image registration
mov=loadtiff(fullfile(home,'raw_data.tif'));
[nrows, ncols, nframes] = size(mov);
movReg=NoRMCorre2(mov); % get registered movie
clear mov
saveastiff(movReg,fullfile(home,'movReg.tif')); % save registered movie
clear movReg

% extract motion traces into MAT file
reg_shifts = returnShifts(home);
save(fullfile(home,'reg_shifts.mat'),'reg_shifts');

%% denoising parameters
mov_in = "movReg.tif";
detr_spacing = 5000;
row_blocks = 4;
col_blocks = 2;
stim_dir = []; % directory with optogenetic stimulation pattern
    
trunc_start = 1; % frame to start denoising
trunc_length = 5000; % length of movie segment to denoise on

%% denoising

run_command = sprintf("source setup.sh\n sbatch denoise.run ""%s"" ""%s"" ""%s"" %d %d %d %d %d ""%s""",...
    home, mov_in, output, detr_spacing, row_blocks, col_blocks,...
    trunc_start-1, trunc_length, stim_dir);

system(run_command);

%% motion correction
moco_command = sprintf("sbatch motion_correction.run ""%s"" ""%s""",...
    home,output);

system(moco_command);

%% blood removal

noise_im = loadtiff(fullfile(output,'Sn_image.tif'));
[ysize, xsize] = size(noise_im);

tStart = tic;
fid = fopen(fullfile(output,'motion_corrected.bin'));                  % open file
tmp = fread(fid, '*float32', 'l');       % uint16, little endian
fclose(fid);                            % close file
L = length(tmp)/(ysize*xsize);
out4 = reshape(tmp, [ysize xsize L]);
clear('tmp');
display(sprintf('Movie loaded successfully. Elapsed time : %.3f s.', toc(tStart)));

smoothing=100; %for high pass filter
movHP = out4 - imfilter(out4, ones(1,1,smoothing)/smoothing, 'replicate');

flucImg = mean(movHP(:,:,1:end-1).*movHP(:,:,2:end), 3);
flucImgS = imfilter(flucImg, fspecial('gaussian', [5 5], 2), 'replicate');

refimg = (flucImgS).^(0.25);

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

bloodmask = uint8(inpoly==0);
mov = out4.*repmat(inpoly==0, [1, 1, nframes]);
saveastiff(bloodmask,fullfile(output,'bloodmask.tif'));

%% background selection
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
