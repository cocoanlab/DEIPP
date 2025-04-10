function question_all = call_question

question_xlsx = fullfile(fileparts(mfilename('fullpath')), 'DEIPP_Questionnaires.xlsx');
question_semantic_list = readtable(question_xlsx, 'Sheet', 'Semantic');
question_therapeutic_list = readtable(question_xlsx, 'Sheet', 'Therapeutic');
question_all.semantic = question_semantic_list.Question;
question_all.therapeutic = question_therapeutic_list.Question;
question_all.semantic = reshape(question_all.semantic, 2, numel(question_all.semantic)/2).';
question_all.therapeutic = reshape(question_all.therapeutic, 2, numel(question_all.therapeutic)/2).';

end

