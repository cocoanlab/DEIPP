function [bootci, bootP] = boot2cip(bootdat, cimethod, bootnull, alpha)

if nargin < 2; cimethod = 'perc'; end
if nargin < 3; bootnull = 0; end
if nargin < 4; alpha = 0.05; end
if iscolumn(bootdat); bootdat = bootdat.'; end
nboot = size(bootdat, 2);

switch cimethod
    case 'norm'
        bootm = mean(bootdat, 2);
        bootstd = std(bootdat, [], 2);
        bootstd(bootstd == 0) = Inf;
        bootci = bootm + norminv([alpha/2 1-alpha/2]) .* bootstd;
        bootZ = (bootm - bootnull) ./ bootstd;
        bootP = 2 * (1 - normcdf(abs(bootZ)));
    case 'perc'
        bootdat_sort = sort(bootdat, 2, 'ascend');
        bootci = bootdat_sort(:, round([nboot*alpha/2 nboot*(1-alpha/2)]));
        bootP = min( ...
            (sum(bootdat_sort <= bootnull, 2)+1) / (nboot+1), ... % Left tail
            (sum(bootdat_sort >= bootnull, 2)+1) / (nboot+1)) .* 2; % Right tail
    otherwise
        error('Unknown method to calculate bootstrap CI and p-value!');
end

end