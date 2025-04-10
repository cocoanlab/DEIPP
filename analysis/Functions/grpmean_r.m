function r_new = grpmean_r(r_orig, grpidx, nodiag)

if nargin < 3; nodiag = true; end

r_new = r_orig;
if nodiag; r_new(1:length(r_new)+1:end) = NaN; end

r_new = splitapply(@(x) mean(x,2,'omitnan'), r_new, findgroups(grpidx(:)).');
r_new = splitapply(@(x) mean(x,1,'omitnan'), r_new, findgroups(grpidx(:)));

end