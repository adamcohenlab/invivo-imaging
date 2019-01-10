clc;clearvars;close all

addpath('quivvia/lib');
addpath('../../../Computer Code/Image Processing');
addpath('../../../Computer Code/Image Processing/NoRMCorre-master')
labpath = labPath();
addpath(fullfile(labpath, 'Labmembers', 'Yoav Adam', 'Scripts', 'NoRMCorre-master'));

data_dir = 'data/';

denoised = double(vm(fullfile(data_dirs{idx},'denoised_15s.tif')));
denoised = permute(denoised,[3 1 2]);
figure(881); clf; moviefixsc(denoised);
refimg = max(denoised(:,:,1000:2000),[],3);

nframes = size(denoised, 3);

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

mov = denoised.*repmat(inpoly==0, [1, 1, nframes]);
bloodmask = uint8(mean(mov,3) ~= 0);
saveastiff(bloodmask,fullfile(data_dirs{idx},'bloodmask.tif'))

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
saveastiff(ff,fullfile(data_dirs{idx},'ff.tif'))
saveastiff(fb,fullfile(data_dirs{idx},'fb.tif'))

exit;