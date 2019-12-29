classdef PCAobj < handle
    %PCAobj defines a PCA algorithm output as a class, as
    %simply as possible, just so it inherits from handle, and has well
    %defined properties which reflect the type of data/analysis that it
    %stores.
    
    properties
        parent = [];
        set = [];
        timestamp = [];
        uvm = [];
        d = [];
        v = [];
        npc = [];
        stats = [];
        more = []; %other outputs which may change during testing
        ICAobj = [];

    end
    
    methods
        function s = saveobj(obj) %function to save without running out of memory
            fnames = fieldnames(obj);
            for i = 1:length(fnames)
                fname = fnames{i};
                if strcmp(fname,'parent')%
                    s.(fname) = [];
                else
                    s.(fname) = obj.(fname);
                end
            end
        end
    end
    
            methods (Static)
        function obj = loadobj(s) %load pcaobj and restore references correctly
            newObj = PCAobj();
            fnames = fieldnames(newObj);
            
            for i = 1:length(fnames)
                fname =fnames{i};
                newObj.(fname)=s.(fname);
            end
            
            if ~isempty(newObj.ICAobj)
                newObj.ICAobj.parent = newObj;
            end
            
                        
            obj = newObj;
        end
    end
end
