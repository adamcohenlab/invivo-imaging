function [obj] = do_binning(obj,b)
%PCA_NORM perform PCA on vm movie object
%   Detailed explanation goes here
    
    if b==1
        return
    else
        Bin = b;
        nrow = obj.rows;
        ncol = obj.cols;
        % chop to make sure movie dimensions divisible by b;
        
        obj = obj((1+mod(nrow,b)):end,:,:);
        obj = obj(:,(1+mod(ncol,b)):end,:);
        [nrow, ncol, nframe2] = size(obj);
        
        mov = obj.data;
        tmp=reshape(mov,Bin,nrow/Bin,Bin,ncol/Bin,nframe2);
        movB=squeeze(mean(mean(tmp,1),3));clear tmp
        obj = vm(movB);
    end


end