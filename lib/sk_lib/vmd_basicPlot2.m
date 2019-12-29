function [h] = vmd_basicPlot2(v,varargin)
%vmd_basicPlot2
%   Plot vmd, with spectrogram and text
%   Simon kheifets 6/29/2008
% add...
%   PSD and autocorrelation of mean trace...
%   statistics abt image...
h.parent = v;

p = inputParser;
p.KeepUnmatched = 1;
addParameter(p,'topdf',0);
addParameter(p,'fixylims',1); %adjust the y limits for mean nad variancce
parse(p,varargin{:});
v2struct(p.Results);
if topdf
    h.fh = figure('Units','inches',...
        'Position',[17 1.5 8.5 11],...
        'Color','w',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperSize',[8.5 11]);
else
    h.fh = figure('Position',[403   100   837   651],'Color','w');
end

ht = 0.1; %height of title
wot = 0.3; %width of object text
hot = 0.5; %height of object text
hims = 0.5*hot; %height of images
htt = hot-hims; %height of the time traces
hsp = 0.4; %height of spectrogram
wsp = 0.7; %width of spectrogram
%axes for the title text
titlesize = 0.15;
htm = 1-ht;
postitle = [0,htm,1,ht];
%axes for the object text
wwot = 0.8*wot;
hhot = 0.8*hot;
wotm = 1-wot;
hotm = 1-hot;
ll = wotm+0.5*(wot-wwot);
bb = 1-hot-ht+0.5*(hot-hhot);
posobjtext = [ll bb wwot hhot];
%axes for images
hh = hims/2;
hhi = 0.8*hh;
ww = wotm;
wwi = 0.8*wotm;
ll =(ww-wwi)/2;
bb1 = 1-ht-hh+0.5*(hh-hhi);
bb2 = 1-ht-hims+0.5*(hh-hhi);
pos1 = [ll,bb1,wwi,hhi];
pos2 = [ll,bb2,wwi,hhi];
%axes for time traces
hh = htt/2;
hhi = 0.8*hh;
ww = wotm;
wwi = 0.8*wotm;
ll =(ww-wwi)/2;
bb1 = 1-hims-ht-hh+0.5*(hh-hhi);
bb2 = 1-hims-ht-htt+0.5*(hh-hhi);
pos3 = [ll,bb1,wwi,hhi];
pos4 = [ll,bb2,wwi,hhi];
%axes for spectrogram and psds

hh = hsp;
ww = wsp;
hhi = 0.8*hsp;
wwi = 0.8*wsp;
ll = 0.5*(ww-wwi);
bb = 0.5*(hh-hhi);
pos_sg = [ll bb wwi hhi];

hh = hsp/2;
ww = 1-wsp;
hhi = 0.8*hh;
wwi = 0.8*ww;
ll = wsp+0.5*(ww-wwi);
bb1 = hh+0.5*(hh-hhi);
bb2= 0.5*(hh-hhi);
pos_ps1 = [ll bb1 wwi hhi];
pos_ps2 = [ll bb2 wwi hhi];

%%title
axcount = 1;
titletext = ['Raw voltage movie summary:' newline v.label];
h.ah(axcount) = subplot('Position',postitle);
%set(gca,'box','on');
h.th = text(0.5, 0.5, titletext, 'Units', 'Normalized',...
    'HorizontalAlignment','Center',...
    'FontUnits','Normalized','FontSize',0.3,...
    'Interpreter','none');
if topdf
    set(h.th,'FontUnits','inches','FontSize',titlesize);
end
set(gca,'Visible','off')

%% info
axcount = axcount+1;
h.ah(axcount)=subplot('Position',posobjtext);
infostring = vmd_toString(v);
h.th2 = text(-0.5,1,infostring,...
    'Interpreter','none',...
    'VerticalAlignment','top',...
    'FontName','FixedWidth');
set(gca,'Visible','off')

%% mean image
axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos1);
imagesc(v.meanimg);
daspect([1 1 1]);
colorbar;
colormap(jet);
title('Mean of each pixel');

%% var image
axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos2);
imagesc(log10(v.varimg));
daspect([1 1 1]);
colorbar;
colormap(jet);
title('Var of each pixel (log scale)');

%% mean trace
axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos3);
plot(v.tvec,v.meantrace,'LineWidth',0.05);
title('Mean of each frame');
xlabel('Time (s)');
ylabel('mean counts');
if fixylims
    ylim([prctile(v.meantrace,.1) max(v.meantrace)]);
end
    
%% var trace
axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos4);
plot(v.tvec,v.vartrace,'LineWidth',0.05);
title('Variance of each frame');
xlabel('time (s)');
ylabel('var (counts^2)');
if fixylims
    ylim([prctile(v.vartrace,.1) max(v.vartrace)]);
end

%% spectrogram
if ~isempty(v.sgrm)
    sgrm = v.sgrm;
    Stot = sgrm.Stot;
    fvecS=sgrm.fvec;
    tvecS = sgrm.tvec;
    fvecSLog=sgrm.fvecLog;
    StotLog = sgrm.StotLog;
    f = sgrm.fvm;
    Svm = sgrm.Svm;


    axcount = axcount+1;
    h.ah(axcount) = subplot('Position',pos_sg);
    %imagesc(log10(Stot),'XData',tvecS,'YData',fvecS);
    cdat = 10*log10(Stot);
    imagesc(tvecS,fvecS,cdat);
    %imagesc(tvecS,fvecSLog, 10*log10(StotLog));
    %colorbar();
    set(gca,'CLim',[prctile(cdat(:),0.5) prctile(cdat(:),99.5)]);
    title('Spectrogram');
    xlabel('time(s)');
    ylabel('f(Hz)');

    %% PSD
    axcount = axcount+1;
    h.ah(axcount) = subplot('Position',pos_ps1);
    [h.xbin, h.ybin] = lowpassplot(f,Svm,300,'bintype','log','Axes',h.ah(axcount));
    xlim([.1 500]);
    ylabel('PSD (counts^2/Hz)');
    title('power density')

    %% PSD Cumulative
    axcount = axcount+1;
    h.ah(axcount) = subplot('Position',pos_ps2);
    semilogx(f(f>10),cumsum(Svm(f>10)));
    xlim([.1 500]);
    title('cumulative power');
    xlabel('frequency (Hz)');
    ylabel('variance (counts^2)');
end
end

