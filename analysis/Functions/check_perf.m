function check_perf(sj_num, opt)

antype = 'cocoan-analysis';
bidsdir = fileparts(fileparts(fileparts(mfilename('fullpath')))); % mfilename: bidsdir/code/~.m
andir = fullfile(bidsdir, 'derivatives', antype);
tbl = tdfread(fullfile(bidsdir, 'participants.tsv'));
sj_id = tbl.participant_id(sj_num,:);
preddir = fullfile(andir, 'predict_rating_func', sj_id);

mdlname = fullfile(preddir, sprintf('%s_predict_rating_modeling_%s_%s_%s%d_tr%d_hout%s_alg%s.mat', ...
    sj_id, opt.dfc, opt.parc, opt.bintype, opt.nbin, opt.ntrain, opt.hout, opt.algorithm));
load(mdlname, 'out');

rcols = call_colors;
idx_train = setdiff(1:size(out.Y,2), out.idx_hout);
lcols = rcols.bipolar(round(linspace(1, size(rcols.bipolar,1), numel(idx_train))), :);

tiledlayout(4, 2);

nexttile(1);
hold on;
for k = 1:numel(idx_train)
    r = cellfun(@(a) a.('r_r'), {out.res.icv(:,idx_train(k)).run}, 'un', false);
    r = cellfun(@(a) mean(a(:, setdiff(1:size(out.Y,2), out.idx_hout)), 1:2, 'omitnan'), r);
    plot(r, 'Color', lcols(k,:));
    [~, i] = max(r);
    plot(i, r(i), '*', 'Color', lcols(k,:));
end
xlabel('Parameter'); ylabel('Run Corr'); box off;
nexttile(3);
hold on;
for k = 1:numel(idx_train)
    r = cellfun(@(a) a.('r_r'), {out.res.icv(:,idx_train(k)).ses}, 'un', false);
    r = cellfun(@(a) mean(a(setdiff(1:size(out.Y,2), out.idx_hout)), 'omitnan'), r);
    plot(r, 'Color', lcols(k,:));
    [~, i] = max(r);
    plot(i, r(i), '*', 'Color', lcols(k,:));
end
xlabel('Parameter'); ylabel('Ses Corr'); box off;
nexttile(5);
hold on;
for k = 1:numel(idx_train)
    r = cellfun(@(a) a.('r_r'), {out.res.icv(:,idx_train(k)).train});
    plot(r, 'Color', lcols(k,:));
    [~, i] = max(r);
    plot(i, r(i), '*', 'Color', lcols(k,:));
end
xlabel('Parameter'); ylabel('All Corr'); box off;
nexttile(7);
hold on;
for k = 1:numel(idx_train)
    r = cellfun(@(a) a.('mse'), {out.res.icv(:,idx_train(k)).train});
    plot(r, 'Color', lcols(k,:));
    [~, i] = min(r);
    plot(i, r(i), '*', 'Color', lcols(k,:));
end
xlabel('Parameter'); ylabel('All MSE'); box off;

nexttile(2);
r = cellfun(@(a) a.('r_r'), {out.res.ocv.run}, 'un', false);
r = cellfun(@(a) mean(a(:, setdiff(1:size(out.Y,2), out.idx_hout)), 1:2, 'omitnan'), r);
hold on;
plot(r, 'Color', 'k');
[~, i] = max(r);
plot(i, r(i), '*', 'Color', 'k');
xlabel('Parameter'); ylabel('Run Corr'); box off;
nexttile(4);
r = cellfun(@(a) a.('r_r'), {out.res.ocv.ses}, 'un', false);
r = cellfun(@(a) mean(a(setdiff(1:size(out.Y,2), out.idx_hout)), 'omitnan'), r);
hold on;
plot(r, 'Color', 'k');
[~, i] = max(r);
plot(i, r(i), '*', 'Color', 'k');
xlabel('Parameter'); ylabel('Ses Corr'); box off;
nexttile(6);
r = cellfun(@(a) a.('r_r'), {out.res.ocv.train});
hold on;
plot(r, 'Color', 'k');
[~, i] = max(r);
plot(i, r(i), '*', 'Color', 'k');
xlabel('Parameter'); ylabel('All Corr'); box off;
nexttile(8);
r = cellfun(@(a) a.('mse'), {out.res.ocv.train});
hold on;
plot(r, 'Color', 'k');
[~, i] = min(r);
plot(i, r(i), '*', 'Color', 'k');
xlabel('Parameter'); ylabel('All MSE'); box off;

set(gcf, 'Position', [361          93        1229         892]);
set(gcf, 'color', 'w');

end