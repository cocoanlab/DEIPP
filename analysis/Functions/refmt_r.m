function r_new = refmt_r(r_orig, nodiag)

if nargin < 2; nodiag = true; end
nodiag = double(logical(nodiag));
yesdiag = double(~nodiag);

if ~isvector(r_orig) && size(r_orig,1)==size(r_orig,2)
    
    % matrix to vector
    r_new = r_orig(triu(true(size(r_orig)), nodiag));
    
elseif isvector(r_orig)
    
    % vector to matrix
    n = 1/2 + sqrt(numel(r_orig)*2 + 1/4) - yesdiag; % if diagonal is included, n = n - 1
    r_new = zeros(n, n);
    r_new(triu(true(size(r_new)), nodiag)) = r_orig;
    r_new = r_new + triu(r_new, 1).';
    
end

end