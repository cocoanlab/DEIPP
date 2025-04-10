%% Save run order.

addpath(genpath(fileparts(mfilename('fullpath'))));
runs = {'REST', 'RATE1', 'RATE2', 'RATE3', 'PREP', 'SPEAK1', 'SPEAK2', 'LISTEN1', 'LISTEN2', 'RELISTEN1', 'RELISTEN2'};
durs = [10, 10, 10, 10, 7, 5, 5, 5, 5, 5, 5] .* 60; % secs
n_subj = 10;
n_visit = 30;
n_run = numel(runs);
marker_mat = false(n_subj, n_visit, n_run);
marker_mat(:, 1, 1) = true;

question_all = call_question;
n_question = size(question_all.semantic, 1);
qmarker_mat = false(n_subj, n_question);
qmarker_mat(:, 1) = true;

save('DEIPP_run_data.mat', 'runs', 'durs', 'marker_mat', 'qmarker_mat');
