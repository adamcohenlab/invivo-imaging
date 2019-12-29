function [h] = fvm_quickPlot(f,varargin)
%fvm_quickplot
%   Plot useful information about filtered movie
%Simon kheifets 6/3/2018
h.parent = f;

p = inputParser;
p.KeepUnmatched = 1;
addParameter(p,'topdf',0);

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
%% define axes' positions
titlesize=0.15; %font size for title
if topdf
    txh0 = 0.1;
    divh0 = 0.4; %how much vert space for the three time traces
    divw0 = 0.6; %how much horizontal space for the four figures
    postitle = [0,1-txh0,1,txh0];
    
    %position for 3 time traces
   
    w0 = 1;
    h0 = (1-txh0)*divh0/3;
    deltaw = 0.9;
    mdeltaw = 1-deltaw;
    deltah = 0.7;
    mdeltah = 1-deltah;
    ww = deltaw*w0;
    hh = deltah*h0;
    for i = 1:3
        lx=mdeltaw*w0/2;
        bx= (i-1)*(h0)+mdeltah*h0/2;
        pos_t{i}=[lx bx ww hh]; %time trace positions
    end
    
    %position for %4 images
    
    nrows = 4;
    ncols = 1;
    nplots = nrows*ncols;
    h0 = (1-txh0)*(1-divh0)/nrows;
    w0 = divw0/ncols;
    deltaw = 0.9;
    mdeltaw = 1-deltaw;
    deltah = 0.8;
    mdeltah = 1-deltah;
    hh = deltaw*h0;
    ww = deltah*w0;
    for i = 1:nrows
        b = 1-txh0-i*h0+mdeltah*h0/2;
        for j=1:ncols
            l = (j-1)*w0+mdeltaw*w0/2;
            pos_f{i,j} = [l b ww hh];
        end
    end
    
    pos{1,1}=pos_f{1};
    pos{1,2}=pos_f{2};
    pos{2,2}=pos_f{3};
    pos{3,2}=pos_f{4};
    pos(:,3)=pos_t;
    
    %position for info text
    ww = 1-divw0;
    hh = 1-divh0-txh0;
    w0 = 0.8*ww;
    h0 = 0.8*hh;
    ll = 1-ww+0.5*(ww-w0);
    bb = 1-txh0-hh+0.5*(hh-h0);
    posobjtext = [ll bb w0 h0];
    
    
else
%     nplots = 9;
    nrows = 3;
    ncols = 3;
    postitle = [0,0.9,1,0.1];
    w0 = 1/ncols;
    h0 = 0.9/nrows;
    ww = (1/ncols)*0.8;
    hh = (0.9/nrows)*0.65;
    for i = 1:nrows
        b = 0.9-i*h0+0.5*(h0-hh);
        for j=1:ncols
            l = (j-1)*w0+0.5*(w0-ww);
            pos{i,j} = [l b ww hh];
        end
    end
end

        

titletext = ['Filtered voltage movie summary: ' newline f.parent.label];
%% title
axcount = 1;
h.ah(axcount) = subplot('Position',postitle);

h.th = text(0.5, 0.5, titletext, 'Units', 'Normalized',...
    'HorizontalAlignment','Center',...
    'FontUnits','Normalized','FontSize',0.3,'Interpreter','none');
set(gca,'Visible','off')

if topdf
    set(h.th,'FontUnits','inches','FontSize',titlesize);
end

%% mean img
axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos{1,1});
imagesc(f.meanimg);
daspect([1 1 1]);
colorbar;
colormap(jet);
title('mean img');

%% var img
axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos{1,2});
imagesc(log10(f.varimg));
daspect([1 1 1]);
colorbar;
colormap(jet);
title('var img');

%% counts
axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos{1,3});
plot(f.tvec,f.meantrace);
title('frame avg');
%xlabel('Time (s)');
ylabel('counts');

%%
% axcount = axcount+1;
% h.ah(axcount) = subplot('Position',pos{2,1});
%imagesc(f.meanimgHP);
% daspect([1 1 1]);
% colorbar;
% colormap(jet);
% title(f.label);

%% hp var img
axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos{2,2});
imagesc(log10(f.varimgHP));
daspect([1 1 1]);
colorbar;
colormap(jet);
title('HP: var img)');

%% hp avg img
axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos{2,3});
plot(f.tvec,f.meantraceHP);
title('HP: frame avg');
%xlabel('Time (s)');
ylabel('counts');
%%
% axcount = axcount+1;
% h.ah(axcount) = subplot('Position',pos{3,1});
% imagesc(f.meanimgLP);
% daspect([1 1 1]);
% colorbar;
% colormap(jet);
% title(f.label);

%% lp var img
axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos{3,2});
imagesc(log10(f.varimgLP));
daspect([1 1 1]);
colorbar;
colormap(jet);
title('LP: var img');

%% lp avg img
axcount = axcount+1;
h.ah(axcount) = subplot('Position',pos{3,3});
plot(f.tvec,f.meantraceLP);
title('LP: frame avg');
xlabel('Time (s)');
ylabel('counts');

%% info text
axcount = axcount+1;
h.ah(axcount)=subplot('Position',posobjtext);
infostring = vmd_toString(f);
h.th2 = text(-0.5,1,infostring,...
    'Interpreter','none',...
    'VerticalAlignment','top',...
    'FontName','FixedWidth');
set(gca,'Visible','off')
end

