function [objout,xShifts,yShifts] = fixMotionEdges(objin,dirpath,varargin)
%FIXMOTIONEDGES load the motion correction data and correct the file
%   objin is a vm object
%   dirpath is the path to the directory which contains
%   'CorrectionMetrics.fig', where the motion information is saved 
%   Simon Kheifets 5/20/2018

 try
     fig = openfig(fullfile(dirpath,'CorrectionMetrics.fig'));
     lxhs = get(get(fig,'Children'),'Children');
     if length(lxhs)==6
        ydata = get(lxhs{6},'YData');
        xShifts = ydata{1};
        yShifts = ydata{2};
     else
         ydata = get(lxhs{2},'YData');
         xShifts = ydata{1};
         ydata = get(lxhs{1},'YData');
         yShifts = ydata{1};
     end
 catch
    fig = openfig(fullfile(dirpath,'RigidCorrectionMetrics.fig'));
    lxhs = get(get(fig,'Children'),'Children');
    ydata = get(lxhs{6},'YData');
    xShifts = ydata{1};
    yShifts = ydata{2};
 end

close(fig);
if nargin==1
    xShifts = xShifts(varargin{1});
    yShifts = yShifts(varargin{1});
end
delta = 0.5;
leftEdge = ceil(prctile(xShifts, 100-delta));
rightEdge = ceil(abs(prctile(xShifts, delta)));
botEdge = ceil(abs(prctile(yShifts, delta)));
topEdge = ceil(prctile(yShifts, 100-delta));
objout = objin((1+topEdge):(end-botEdge), (1+leftEdge):(end-rightEdge),:);

end

