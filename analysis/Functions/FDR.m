function [pID,pN] = FDR(p,q)
% :Usage:
% ::
%
%     pt = FDR(p,q)
%
% :Inputs:
%
%   **p:**
%        vector of p-values
%
%   **q:**
%        False Discovery Rate level
%
% :Outputs:
%
%   **pID:**
%        p-value threshold based on independence or positive dependence
%
%   **pN:**
%        Nonparametric p-value threshold
%
% ..
%    % @(#)FDR.m	1.3 Tom Nichols 02/01/18
%
%    The checking code below was added by Tor Wager
% ..

    p(isnan(p)) = [];

    p = sort(p(:));
    V = length(p);
    I = (1:V)';

    cVID = 1;
    cVN = sum(1./(1:V));

    pID = p(max(find(p<=I/V*q/cVID)));
    pN = p(max(find(p<=I/V*q/cVN)));

    return


