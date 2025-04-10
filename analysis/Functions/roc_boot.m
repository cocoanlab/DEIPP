function [ci, names] = roc_boot(input_vals, binary_outcome, thr, class_thr, verbose)
% :Usage:
% ::
%
%     [ci, names] = roc_boot(input_vals, binary_outcome, thr, [verbose flag])
%
% Returns bootstrapped 95% confidence intervals for sensitivity,
% specificity, and PPV at a given threshold.
%
% :Examples:
% ::
%
%    thr = ROC.class_threshold
%    input_vals = input(ind);
%    binary_outcome = outcome(ind);


if nargin < 4, verbose = 1; end

bootfun = @(input_vals, binary_outcome) classification_stats(input_vals, binary_outcome, thr, class_thr);

sens_spec_ppv_auc = bootstrp(10000, bootfun, input_vals, binary_outcome);

% 95% CIs
for i = 1:size(sens_spec_ppv_auc,2)
    ci{i} = [prctile(sens_spec_ppv_auc(:, i), 2.5)  prctile(sens_spec_ppv_auc(:, i), 97.5)];
end

names = {'95% CI for sensitivity' '95% CI for specificity' '95% CI for PPV', '95% CI for AUC'};

if verbose
    ssp = bootfun(input_vals, binary_outcome);
    
    for i = 1:size(sens_spec_ppv_auc,2)
        fprintf('%s\t%.0f%% +- (%.0f-%.0f%%)\n', names{i}, ssp(i)*100, ci{i}(1)*100, ci{i}(2)*100);
    end
    
end


end % main function


function sens_spec_ppv_auc = classification_stats(input_vals, binary_outcome, thr, class_thr)

% sens_spec_ppv
wh = input_vals >= class_thr;

tpr = sum(wh(binary_outcome)) ./ sum(binary_outcome);
fpr = sum(wh(~binary_outcome)) ./ sum(~binary_outcome);

npos = sum(tpr .* binary_outcome);
nfp = sum(fpr .* ~binary_outcome);
ppv = npos ./ (npos + nfp);

sens_spec_ppv_auc = [tpr 1-fpr ppv];

% AUC
[tpr, fpr] = deal(zeros(size(thr)));

indx = 1;
for x = thr
    wh = input_vals >= x;
    
    tpr(indx) = sum(wh(binary_outcome)) ./ sum(binary_outcome);
    fpr(indx) = sum(wh(~binary_outcome)) ./ sum(~binary_outcome);
    
    indx = indx + 1;
end

auc = calc_auc(fpr, tpr);

sens_spec_ppv_auc = [sens_spec_ppv_auc auc];

end

function auc = calc_auc(fpr, tpr)

[u, wh] = unique(fpr);
u2 = tpr(wh);

% fix for AUC = 1 if no overlap; triangle method not perfectly accurate
% here.
if any(u == 0 & u2 == 1), auc = 1; return, end
if isequal(u, 1); auc = 0; return; end

for i = 2:length(u)
    
    xdiff = u(i) - u(i - 1);
    ydiff = u2(i) - u2(i - 1);
    a(i) = xdiff * u2(i - 1) + xdiff * ydiff / 2;  % area of rect + area of triangle
    
end


auc = sum(a);

end
