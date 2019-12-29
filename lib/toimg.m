function out = toimg(mat, nr, nc);

% function out = toimg(mat, nr, nc);
% converts a matrix of size [nr*nc, nt] to a movie of size [nr, nc, nt];
%
% modified from code by Vicente Parot

out = reshape(mat, nr, nc, []);