function [out, scaleImgs] = SeeResiduals(mov, refSig, remOffset); 

% function [out, scaleImgs] = SeeResiduals(in, refSig, remOffset); 
% project out the component of a movie with timetrace given by refSig.
% size(Movie) = [Y, X, T];
% size(refSig) = [N, T], where N is the number of reference signals to
% remove.
% if refSig == 0, then only the average image is subtracted from each
% frame.
% 
% remOffset: Boolean variable.  If 1, add a row of 1's to the refSig, to
% remove constant offsets at each pixel.  If 0, do not.  If omitted,
% default = 1.
%
% out = corrected movie
% scaleImgs = images of the coefficients for each temporal trace removed.
% size(scaleImages) = [X, Y, N], or [X, Y, N+1] if remOffset == 1.
%
% AEC 13 May 2015

if isa(mov,'vm') % 2017 Sami, VP, added vm compatibility
    mov = mov.toimg.data;
end

if nargin == 2
    remOffset = 1;  % Default: if not specified whether to remove offset, then subtract off an average image.
end
if ~isscalar(remOffset)
    'remOffset must be a scalar (zero or nonzero)'
    return
end
% 
[ySize, xSize, nFrame] = size(mov);

if refSig == 0 % Then only subtract the average image from each frame.
    avgImg = mean(mov, 3);
    out = mov - repmat(avgImg, [1 1 nFrame]);
    scaleImgs = avgImg;
else  %Project out the time-trace(s) in refSig.
    % Make sure refSig is the right size and orientation
    [refY, refX] = size(refSig);
    if refY == nFrame && refX ~= nFrame
        refSig = refSig';
    elseif refX ~= nFrame
        'Size of refSig must be [N, # of frames]'
        return
    end
    
    if remOffset
        dRef = bsxfun(@minus, refSig, mean(refSig, 2));  % remove the DC offset from refSig because we will account for offset separately.
        I = [ones(1, nFrame); dRef];  % Add a row of ones to measure the DC offset at each pixel.
    else
        I = refSig;
    end

    movVec = tovec(mov);
    C = movVec*I'*inv(I*I');  % Linear algebra to find the pseudo-inverse
    scaleImgs = toimg(C, ySize, xSize);  % convert weights to images.

    residMovVec = movVec - C*I;  % Look at the residuals
    out = toimg(residMovVec, ySize, xSize);
end;


