function [handles] = stackImshow(obj,varargin)
%ah = STACKIMSHOQ(vm,'parameter name',parameter value) plot a grid of movie frames
%   Detailed explanation goes here
% obj: vm movie object
% frames: indices of frames (if [], then does all frames)
% maintitle: string for main title;
% subtitles: cell of strings for subtitles
% hint: use secret function sprintfc('N=%i',1:30) to quickly generate
% subtitles
%
% Simon Kheifets April 2018

%titletext = 'title text';


p = inputParser;
p.KeepUnmatched = 1;%don't give error when unused parameters appear
addParameter(p,'frames',[]);    %indices of frames (if [], then does all frames)
addParameter(p,'maintitle',[]); %string for main title;
addParameter(p,'subtitles',[]); %cell of strings for subtitles
addParameter(p,'hframe',[]);    %parent axes if you want to embed it
addParameter(p,'tsp',0.06);     %space for title
addParameter(p,'ncols',[]);     %number of columns
addParameter(p,'addcolorbar',0);%put colorbar
addParameter(p,'rbcmap',1);     %use redwhiteblue colormap

if nargin >1 %
    parse(p,varargin{:});
else parse(p);
end
v2struct(p.Results);

if isempty(frames)
    frames = 1:obj.frames;
end
nFrames = length(frames);
if nFrames>100
    error('too many frames!');
    return
end

if isempty(subtitles)
    dosubtitles = 0;
else
    dosubtitles = 1;
end

if isempty(maintitle)
    domaintitle=0;
    tsp=0;
else
    domaintitle=1;
end

rAspect = obj.cols/obj.rows; %aspect ratio of movie frames
if isempty(ncols)
    nCols = max(1,floor(sqrt(nFrames/rAspect)));
else
    nCols = ncols;
end
nRows = ceil(nFrames/nCols);





%set dimensions
if isempty(hframe)
    parentPos = [0 0 1 1];
else
    parentPos = get(hframe,'Position');
    %set(hframe,'Visisble','off');
end
L=parentPos(1); B=parentPos(2); W=parentPos(3); H=parentPos(4);
w = W/nCols; %width between frames
h = H*(1-tsp)/nRows; %height between frames
gx = W*0.02; %gap between frames (horizontal)
gy = H*0.02; %gap between frames (vertical)

%plot the images
for i = 1:nFrames
    iR = ceil(i/nCols); %which row
    iC = i-(iR-1)*nCols; %which column
    iL = L+(w*(iC-1));
    iB = B+(nRows-iR)*h;

    ah(i) = axes('Position',[iL+gx/2 iB+gy/2 (w-gx) (h-gy)]);
    
    set(gca,'box','on');
    f = frame(obj,frames(i));
    imh = imshow(f,[]);
    
    
    if addcolorbar==1
        colorbar();
    end
    set(gca,'CLim', max(abs(f(:)))*[-1 1]);
    colormap(bluewhitered)
    if dosubtitles
        ahst(i) = axes('Position',[iL+gx/2 iB+h-gy/2 w-gx gy]);
        set(gca,'xtick',[],'ytick',[])
        set(gca,'box','on');
        set(ahst(i),'Visible','off');
        thst(i) = text(0, 0, subtitles{i}, 'Units', 'Normalized',...
            'HorizontalAlignment','Left','VerticalAlignment','Bottom','FontUnits','Normalized','FontSize',1);
        handles.ahsubtitles = ahst;
        handles.thsubtitles = thst;
    end
    
end

if domaintitle
% plot the 'title'
%aht = axes('Position', [L+gx/2 B+H-tsp+gy/2 W-gx H*tsp-gy]);
aht = axes('Position', [gx/2 B+H-tsp+gy/2 1-gx H*tsp-gy]);
set(gca,'box','on');
tht = text(0.5, 0.5, maintitle, 'Units', 'Normalized',...
    'HorizontalAlignment','Center','FontUnits','Normalized','FontSize',1);
set(aht,'Visible','off')
handles.thtitle = tht;
end

handles.ahimages = ah;





