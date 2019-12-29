function [fovpath] = fovPath(uprDatas)
%UNTITLED5 Summary of this function goes here
% example of uprdata array:
%  uprData = [...
%     48,2018,05,15,4,1,3,...
%     51,2018,05,31,1,1,2,...
%     51,2018,05,31,1,1,3,...
%     52,2018,05,31,1,1,1,...
%     52,2018,05,31,1,1,2]; %ivq# ,y,m,d, session,slice,fov,
for i = 1:size(uprDatas,1)
    uprData = uprDatas(i,:);
    s1 = sprintf('IVQ%i',uprData(1));
    s2 = sprintf('%i-%0.2i-%0.2i_IVQ%i-S%i',uprData([2 3 4 1 5]));
    s3 = sprintf('slice%i',uprData(6));
    s4 = sprintf('FOV%i',uprData(7));
    fovpath{i} = fullfile(s1,s2,s3,s4);
    if size(uprDatas,1)==1
        fovpath = fovpath{1};
    end

end

    
end

