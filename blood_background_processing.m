clc;clearvars;close all

addpath('quivvia/lib');
addpath('../../../Computer Code/Image Processing');
addpath('../../../Computer Code/Image Processing/NoRMCorre-master')
labpath = labPath();
addpath(fullfile(labpath, 'Labmembers', 'Yoav Adam', 'Scripts', 'NoRMCorre-master'));

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
regalbasepath = '/n/regal/cohen_lab/mxie/data/';
data_dirs = cellfun(@(x) fullfile(regalbasepath,x,'Final'),fovpaths,'UniformOutput',false);
nMovs = size(data_dirs,2);

parfor_progress(nMovs);
for idx = 1:nMovs
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
    parfor_progress;
end
parfor_progress(0);