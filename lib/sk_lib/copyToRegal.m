function [po status msg] = copyToRegal(ri,ro,rel)
%ri 'root' path in
%ro 'root' path out
%rel struct of paths to be copied relative to ro and ri
np = length(rel); %number of paths
pi = cellfun(@(x) fullfile(ri,x),rel,'UniformOutput',false);
po = cellfun(@(x) fullfile(ro,x),rel,'UniformOutput',false);

pc=parcluster('local');
pc.NumWorkers=10;

parpool('local', 10);
parfor_progress(np);
parfor i = 1:np
    fprintf('copying from: %s\n',pi{i});
    if ~isdir(pi{i})
        error(['could not find' pi{i}]);
        
    end
    if ~isdir(po{i})
        mkdir(po{i});
    end
    [status{i}, msg{i}]=copyfile(pi{i},po{i});
    parfor_progress;
end
parfor_progress(0);