function [ics, mixmat, sepmat] = sorted_ica2(pcs,n)
%SORTED_ICA returns sorted independent components.
%   [ics, mixmat, sepmat] = sorted_ica(pcs,n) computes ICA using
%   fastica(pcs,n). Independent components are then sorted according to the
%   decreasing absolute skewness of their median-highpassed versions. A
%   temporal, 30 point median filtered version of each trace is subtracted
%   for sorting. This removes some artifacts that otherwise survive the
%   skewness sorting. The unfiltered traces are then sorted in the computed
%   order and returned. The mixing matrix and the separation matrix are
%   also returned with the same ordering.
%
%   ics = pcs*sepmat'
%   If n == size(pcs, 2) (i.e. # of ics == # of pcs) then
%   pcs = ics*mixmat'
%
%
%   2014 Vicente Parot
%   Cohen Lab - Harvard University
%
%   Modifications made by SK Apr 2018:
%   -add 'g' parameter
%   -don't transpose input
%   -keep output the same orientation as input
%
    thistic = tic;
    disp 'computing ica ...'
    [ics0,mixmat0,sepmat0] = fastica(pcs,... % ica
        'numofic',n,...
        'verbose','off',...
        'g', 'skew');
    ics0 = ics0';
    % apply sorting on median removed version. sort, keep and return the
    % unfiltered signals.
    icskews = skewness(ics0);
%     icskews = skewness(ics0 - medfilt2(ics0,[30 1],'symmetric'));
    [~,iidx] = sort(abs(icskews),'descend'); % sort median filtered ic skews, get indexing values
    ics = bsxfun(@times,ics0(:,iidx),sign(icskews(iidx))); % sort ics and invert negative skews
    mixmat = bsxfun(@times,mixmat0(:,iidx),sign(icskews(iidx))); % mixmat
    sepmat = bsxfun(@times,sepmat0(iidx,:),sign(icskews(iidx)')); % sepmat
    disp(['ica took ' num2str(toc(thistic)) ' s']);
    ics = ics';
end
