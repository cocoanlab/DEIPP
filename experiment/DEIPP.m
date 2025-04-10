%% SETUP : Basic parameter

clear;
close all;
Screen('CloseAll');

global theWindow lb1 rb1 lb2 rb2 H W scale_W anchor_lms korean alpnum space special bgcolor white orange red;

test_mode = true;
check_recording = true;

if test_mode
    USE_BIOPAC = false;
    show_cursor = false; 
    screen_mode = 'full';
    disp('***** TEST mode *****');
else
    USE_BIOPAC = true;
    show_cursor = false;
    screen_mode = 'full';
    disp('***** EXPERIMENT mode *****');
end

[~, hostname] = system('hostname');
switch strtrim(hostname)
    case 'JJ-macpro.local'
        basedir = '/Users/jaejoonglee/Dropbox/Paperwork/DEIPP_sync/experiment/DEIPP_ex';
        inputdev.HostAudioAPIName = 'Core Audio';
        inputdev.DeviceName = 'JaeJoong Lee�� AirPods'; %'���� ����ũ'
        outputdev.HostAudioAPIName = 'Core Audio';
        outputdev.DeviceName = 'JaeJoong Lee�� AirPods'; %'���� ���'
    case 'Cocoanui-iMac-2.local'
        basedir = '/Users/deipp/Dropbox/Paperwork/DEIPP_sync/experiment/DEIPP_ex';
        inputdev.HostAudioAPIName = 'Core Audio';
        inputdev.DeviceName = '���� ����ũ';
        outputdev.HostAudioAPIName = 'Core Audio';
        outputdev.DeviceName = '���� ���';
    case 'DESKTOP-3TBJ304'
        basedir = 'C:\Users\Cocoanlab_WL01\Dropbox\Paperwork\DEIPP_sync\experiment\DEIPP_ex';
        inputdev.HostAudioAPIName = 'ASIO';
        inputdev.DeviceName = 'ASIO4ALL v2';
        outputdev.HostAudioAPIName = 'ASIO';
        outputdev.DeviceName = 'ASIO4ALL v2';
end

cd(basedir);
addpath(genpath(basedir));

exp_scale = {'cont_int_vas', 'overall_alertness', 'cont_threat_vas'};
main_scale = {'cont_int_vas', 'cont_threat_vas'};


%% SETUP : Check subject info

subjID = input('\nSubject ID? ', 's');
subjID = strtrim(subjID);

subjnum = input('\nSubject number? ');

subjvis = input('\nVisit number? ');

subjrun = input('\nRun number? ');

subjquestion = input('\nQuestion number? ');

if isempty(subjID) || isempty(subjnum) || isempty(subjvis) || isempty(subjrun) || isempty(subjquestion)
    error('Wrong number. Break.')
end


%% SETUP : Load randomized run data and Compare the markers

markerfile = fullfile(basedir, 'DEIPP_run_data.mat');
load(markerfile, 'runs', 'durs', 'marker_mat', 'qmarker_mat');
[visitmarker, runmarker] = find(squeeze(marker_mat(subjnum, :, :)));
questionmarker = find(qmarker_mat(subjnum,:));

if visitmarker ~= subjvis
    cont_or_not = input(['\nThe visit number is inconsistent with the latest progress. Continue?', ...
        '\n1: Yes.  ,   2: No, break.\n:  ']);
    if cont_or_not == 1
        visitmarker = subjvis;
    elseif cont_or_not == 2
        error('Break.')
    else
        error('Wrong number. Break.')
    end
end

if runmarker ~= subjrun
    cont_or_not = input(['\nThe run number is inconsistent with the latest progress. Continue?', ...
        '\n1: Yes.  ,   2: No, break.\n:  ']);
    if cont_or_not == 1
        runmarker = subjrun;
    elseif cont_or_not == 2
        error('Break.')
    else
        error('Wrong number. Break.')
    end
end

if questionmarker ~= subjquestion
    cont_or_not = input(['\nThe Question number is inconsistent with the latest progress. Continue?', ...
        '\n1: Yes.  ,   2: No, break.\n:  ']);
    if cont_or_not == 1
        questionmarker = subjquestion;
    elseif cont_or_not == 2
        error('Break.')
    else
        error('Wrong number. Break.')
    end
end


%% SETUP : Save data in first

savedir = fullfile(basedir, 'Data');

nowtime = clock;
subjtime = sprintf('%.2d%.2d%.2d', nowtime(1), nowtime(2), nowtime(3));

data.subject = subjID;
data.datafile = fullfile(savedir, sprintf('%s_%s_sub_%.3d_visit_%.2d_run_%.2d.mat', subjtime, subjID, subjnum, visitmarker, runmarker));
data.version = 'DEIPP_cocoanlab_20220106';
data.starttime = datestr(clock, 0);
data.starttime_getsecs = GetSecs;

save(data.datafile, 'data');


%% SETUP : Paradigm

S.type = runs{runmarker};
S.dur = durs(runmarker);
if test_mode; S.dur = S.dur ./ 60; end
S.stimtext = '+';
S.changetime = 1;

data.dat.type = S.type;
data.dat.duration = S.dur;
data.dat.stimtext = S.stimtext;
data.dat.exp_scale = exp_scale;
data.dat.main_scale = main_scale;

rating_types = call_ratingtypes;
question_all = call_question;

postrun_start_t = 2; % postrun start waiting time.
postrun_end_t = 2; % postrun questionnaire waiting time.

wh_tx = [10, 20, 21:30];


%% SETUP : BIOPAC

if ismember(S.type, {'PREP', 'RELISTEN1', 'RELISTEN2'})
    USE_BIOPAC = false;
end

if USE_BIOPAC
    channel_n = 1;
    biopac_channel = 0;
    ljHandle = BIOPAC_setup(channel_n);
end


%% SETUP : Display the run order

runs_for_display = runs;
runs_for_display{runmarker} = sprintf('[[%s]]', runs_for_display{runmarker});
fprintf('\n\n');
fprintf('Runs: %s\n\n', string(join(runs_for_display)));
input('To continue, press any key.');


%% SETUP : Screen

PsychDefaultSetup(1);
screens = Screen('Screens');
window_num = screens(end);
if ~show_cursor
    HideCursor;
end

Screen('Preference', 'SkipSyncTests', 1);

[window_width, window_height] = Screen('WindowSize', window_num);
switch screen_mode
    case 'full'
        window_rect = [0 0 window_width window_height]; % full screen
    case 'semifull'
        window_rect = [0 0 window_width-100 window_height-100]; % a little bit distance
    case 'middle'
        window_rect = [0 0 window_width/2 window_height/2];
    case 'small'
        window_rect = [0 0 400 300]; % in the test mode, use a little smaller screen
end

% size
W = window_rect(3); % width
H = window_rect(4); % height

lb1 = W/4; % rating scale left bounds 1/4
rb1 = (3*W)/4; % rating scale right bounds 3/4

lb2 = W/3; % new bound for or not
rb2 = (W*2)/3;

scale_W = (rb1-lb1).*0.1; % Height of the scale (10% of the width)

anchor_lms = [0.014 0.061 0.172 0.354 0.533].*(rb1-lb1)+lb1;

% font
fontsize = 33;
Screen('Preference', 'TextEncodingLocale', 'ko_KR.UTF-8');

% color
bgcolor = 50;
white = 255;
red = [158 1 66];
orange = [255 164 0];

% open window
theWindow = Screen('OpenWindow', window_num, bgcolor, window_rect); % start the screen
Screen('Textfont', theWindow, '-:lang=ko');
Screen('TextSize', theWindow, fontsize);

% get font parameter
[~, ~, wordrect1, ~] = DrawFormattedText(theWindow, double('��'), lb1-30, H/2+scale_W+40, bgcolor);
[~, ~, wordrect2, ~] = DrawFormattedText(theWindow, double('p'), lb1-30, H/2+scale_W+40, bgcolor);
[~, ~, wordrect3, ~] = DrawFormattedText(theWindow, double('p '), lb1-30, H/2+scale_W+40, bgcolor);
[~, ~, wordrect4, ~] = DrawFormattedText(theWindow, double('^'), lb1-30, H/2+scale_W+40, bgcolor);
[korean.x, korean.y, alpnum.x, alpnum.y, space.x, space.y, special.x, special.y] = deal(...
    wordrect1(3)-wordrect1(1), wordrect1(4)-wordrect1(2), ... % x = 36, y = 50
    wordrect2(3)-wordrect2(1), wordrect2(4)-wordrect2(2), ... % x = 25, y = 50
    wordrect3(3)-wordrect3(1) - (wordrect2(3)-wordrect2(1)), wordrect3(4)-wordrect3(2), ... % x = 12, y = 50
    wordrect4(3)-wordrect4(1), wordrect4(4)-wordrect4(2)); % x = 19, y = 50

Screen(theWindow, 'FillRect', bgcolor, window_rect); % Just getting information, and do not show the scale.
Screen('Flip', theWindow);


%% SETUP: Input setting (for Mac and test)

Exp_key = [];
Par_key = [];
Scan_key = [];

if ismember(S.type, {'REST', 'SPEAK1', 'SPEAK2', 'LISTEN1', 'LISTEN2', 'RELISTEN1', 'RELISTEN2'})
    InitializePsychSound;
    padev = PsychPortAudio('GetDevices');
    padev_input = padev([padev.NrInputChannels] > 0);
    padev_input = padev_input(strcmp({padev_input.DeviceName}, inputdev.DeviceName) & strcmp({padev_input.HostAudioAPIName}, inputdev.HostAudioAPIName));
    padev_output = padev([padev.NrOutputChannels] > 0);
    padev_output = padev_output(strcmp({padev_output.DeviceName}, outputdev.DeviceName) & strcmp({padev_output.HostAudioAPIName}, outputdev.HostAudioAPIName));
end


%% SETUP: Introduction, Localizer, ANC, Sound test

if strcmp(S.type, 'REST')

    msgtxt = ['�ȳ��ϼ���, �����ڴ�. ''�� ��� ������ ���� ���� ���� ����''�� �������ּż� �����մϴ�.\n\n', ...
        '�� �Կ��� ������ ���� �� 4���� ������ �����Ǿ� ������, �� 1�ð� 30���� �ҿ�˴ϴ�.'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
    inst_img = imread(fullfile(basedir, 'Pictures', 'Instruction_001.png'));
    inst_img_WH = size(inst_img);
    inst_img_WH_resized = inst_img_WH .* (W*3/4 / inst_img_WH(2));
    inst_img_rect = [W/2 - inst_img_WH_resized(2)/2, H*3/4 - inst_img_WH_resized(1)/2, W/2 + inst_img_WH_resized(2)/2, H*3/4 + inst_img_WH_resized(1)/2];
    Screen('PutImage', theWindow, inst_img, inst_img_rect);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    msgtxt = ['�����ڴ��� �Կ��� �ռ� �ִ��� ���� �ڼ��� ���Ͻ� ��, ������ ���� ������ �������ֽñ� �ٶ��,\n' ...
        '�Ӹ��� �����̰ų� �ῡ ���� �ʵ��� ������ �������ֽñ� �ٶ��ϴ�.'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    msgtxt = ['����, �����ڴ��� �Ӹ� ��ġ �ľ��� ���� ���� �Կ��� �����ϰڽ��ϴ�.\n\n', ...
        '�����ڴ��� ������ ���� ��Ÿ���� ȭ�� �߾��� + ǥ�ø� �����ϸ鼭 ����� ��ø� �˴ϴ�.'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    msgtxt = '+';
    DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    msgtxt = ['�������δ�, �Կ� �� ���� ���� �׽�Ʈ�Դϴ�.\n\n', ...
        '�����ڴ��� ȭ�� �߾��� + ǥ�ø� �����ϸ鼭, �������� �ȳ��� ��ٷ��ֽñ� �ٶ��ϴ�.'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    % ANC calibration -> Start dummy scan (5 min) -> ANC learning
    msgtxt = '+';
    DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    % After applying ANC, check sound volume of experimenter and participant
    msgtxt = '��ĳ���� ������ ����� �پ�������, �����ڿ� �������� �Ҹ� ũ�Ⱑ ������ �� Ȯ���մϴ�.';
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    % Finish testing
    msgtxt = ['�׽�Ʈ�� �Ϸ�Ǿ����ϴ�. �����ڴ� �׽�Ʈ ��ĵ�� �����մϴ�.\n\n', ...
        '���ݺ��� �� ������ ���۵˴ϴ�.'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

elseif strcmp(S.type, 'SPEAK1') && check_recording

    msgtxt = ['�Կ��� �ռ�, �Ҹ� �Է�-��� �׽�Ʈ�� �����ϰڽ��ϴ�.\n\n', ...
        '�����ڴ��� ȭ�� �߾��� + ǥ�ø� �����ϸ鼭, �������� �ȳ��� ��ٷ��ֽñ� �ٶ��ϴ�.'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    % Run scan
    msgtxt = '+';
    DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    % Record conversation between experimenter and participant
    msgtxt = ['�����ڿ� �������� ������ ��ȭ ������ �׽�Ʈ�մϴ�.\n\n', ...
        'ȭ�鿡 + ǥ�ð� ��Ÿ���� �������� ������ ������ּ���!']; % ask name and then age
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    pahandle = PsychPortAudio('Open', padev_input.DeviceIndex, 2, 1, [], padev_input.NrInputChannels);
    data.dat.temp_recording_struct = PsychPortAudio('GetStatus', pahandle);
    data.dat.temp_recording_totaldur = 16;
    data.dat.temp_recording_dummydur = 6;
    PsychPortAudio('GetAudioData', pahandle, data.dat.temp_recording_totaldur);
    Screen('Flip', theWindow); % blank screen
    data.dat.temp_recording_starttime = PsychPortAudio('Start', pahandle, 0, 0, 1);
    start_t = GetSecs;
    while true
        cur_t = GetSecs;
        if cur_t - start_t >= data.dat.temp_recording_dummydur
            break
        end
    end
    msgtxt = '+';
    DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true
        cur_t = GetSecs;
        if cur_t - start_t >= data.dat.temp_recording_totaldur
            break
        end
    end
    temp_rec_data = PsychPortAudio('GetAudioData', pahandle);
    temp_rec_data = temp_rec_data(1,:); % Only Lt channel - less noisy
    data.dat.temp_recording_endtime = GetSecs;
    PsychPortAudio('Stop', pahandle);
    PsychPortAudio('Close');
    data.dat.temp_recording_file = data.datafile;
    data.dat.temp_recording_file = strrep(data.dat.temp_recording_file, savedir, fullfile(basedir, 'Recording'));
    data.dat.temp_recording_file = strrep(data.dat.temp_recording_file, '.mat', '_testrecord.wav');
    psychwavwrite(temp_rec_data.', data.dat.temp_recording_struct.SampleRate, data.dat.temp_recording_file);

    % Do post-processing immediately
    msgtxt = ['���� ���� ���� �� ��ĳ�� ���� ���� ���Դϴ�...\n\n', ...
        '�����ڴ��� ���ݸ� ��ٷ� �ֽñ� �ٶ��ϴ�.'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    % Check sound volume of processed recording data (if needed, adjust the Laptop or Line 1 volume)
    msgtxt = ['�Ϸ�Ǿ����ϴ�. ���� ������ ����մϴ�.\n\n', ...
        '�����ڴ��� ������ ��� �Ҹ� ũ�Ⱑ ������ �� Ȯ���մϴ�.'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end
    
    data.dat.temp_recording_proc_file = strrep(data.dat.temp_recording_file, '_testrecord.wav', '_testrecord_p.wav');
    while true 
        if exist(data.dat.temp_recording_proc_file, 'file')
            break
        else
            pause(1);
        end
    end
    while true
        try
            [temp_rec_proc_data, temp_rec_proc_freq] = psychwavread(data.dat.temp_recording_proc_file);
            if size(temp_rec_proc_data,1) / temp_rec_proc_freq >= data.dat.temp_recording_totaldur - 0.1
                break;
            end
        catch
            pause(1);
        end
    end
    if size(temp_rec_proc_data,2) == 1 && padev_output.NrOutputChannels == 2
        temp_rec_proc_data = repmat(temp_rec_proc_data, 1, 2);
    end
    temp_rec_proc_data = temp_rec_proc_data(ceil(temp_rec_proc_freq * data.dat.temp_recording_dummydur)+1 : end, :); % Remove dummy 6 seconds
    pahandle = PsychPortAudio('Open', padev_output.DeviceIndex, 1, 0, temp_rec_proc_freq, padev_output.NrOutputChannels);
    PsychPortAudio('FillBuffer', pahandle, temp_rec_proc_data.');
    PsychPortAudio('Start', pahandle, 0, 0, 1);
    msgtxt = '+';
    DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end
    PsychPortAudio('Stop', pahandle);
    PsychPortAudio('Close');

    % Finish testing
    msgtxt = ['�׽�Ʈ�� �Ϸ�Ǿ����ϴ�. �����ڴ� �׽�Ʈ ��ĵ�� �����մϴ�.\n\n', ...
        '���ݺ��� �� ������ ���۵˴ϴ�.'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

end


%% MAIN : Explain scale

switch S.type


    case 'REST'

        msgtxt = sprintf(['�� �Կ��� %d ��° ������ ''������ �Կ�'' �Դϴ�.\n\n', ...
            '������ �Կ� ����, �����ڴ��� ������ ���� ��Ÿ���� ȭ�� �߾��� + ǥ�ø� �����ϸ鼭 ����� ��ø� �˴ϴ�.'], runmarker);
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end

        msgtxt = '+';
        DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end

        % Explain scale with visualization
        msgtxt = '�� �Կ� ������ ���� ������, ������ ���� �� �򰡰� �̷�����ϴ�.';
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, orange, [], [], [], 2);
        ratetype = strcmp(rating_types.alltypes, exp_scale{2});
        DrawFormattedText(theWindow, rating_types.prompts{ratetype}, 'center', 300, white, [], [], [], 2);
        [lb, rb, start_center] = draw_scale(exp_scale{2});
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end

        % Initial position
        if start_center
            SetMouse((rb+lb)/2,H/2); % set mouse at the center
        else
            SetMouse(lb,H/2); % set mouse at the left
        end

        % Get ratings
        while true % Button
            msgtxt = '�����ڴ��� ����� �� ����� �����Ͻ� ��, ������ ������ ��ư�� �����ֽñ� �ٶ��ϴ�.';
            DrawFormattedText(theWindow, double(msgtxt), 'center', 100, orange, [], [], [], 2);
            DrawFormattedText(theWindow, rating_types.prompts{ratetype}, 'center', 300, white, [], [], [], 2);

            [lb, rb, start_center] = draw_scale(exp_scale{2});

            [x,~,button] = GetMouse(theWindow);
            if x < lb; x = lb; elseif x > rb; x = rb; end
            if button(1); while button(1); [~,~,button] = GetMouse(theWindow); end; break; end
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end

            Screen('DrawLine', theWindow, orange, x, H/2, x, H/2+scale_W, 6);
            Screen('Flip', theWindow);
        end

        % Freeze the screen 0.5 second with red line
        Screen('DrawLine', theWindow, red, x, H/2, x, H/2+scale_W, 6);
        Screen('Flip', theWindow);

        freeze_t = GetSecs;
        while true
            freeze_cur_t = GetSecs;
            if freeze_cur_t - freeze_t > 0.5
                break
            end
        end


    case 'RATE1'

        msgtxt = sprintf(['�� �Կ��� %d ��° ������ ''���� ��'' �Դϴ�.\n\n', ...
            '�Կ� ���� �����ڴ��� ������ ���� ô�� ���� ���� �������� ������ ���⸦ ���������� �����ֽø� �˴ϴ�.'], runmarker);
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end

        % Explain scale with visualization
        msgtxt = '���� �� ����';
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, orange, [], [], [], 2);
        ratetype = strcmp(rating_types.alltypes, exp_scale{1});
        DrawFormattedText(theWindow, rating_types.prompts{ratetype}, 'center', 300, white, [], [], [], 2);
        [lb, rb, start_center] = draw_scale(exp_scale{1});
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end

        % Initial position
        if start_center
            SetMouse((rb+lb)/2,H/2); % set mouse at the center
        else
            SetMouse(lb,H/2); % set mouse at the left
        end

        % Get ratings
        while true % Space
            msgtxt = '�����ڴ��� ����� �� ����� �����Ͻ� ��, ������ ������ ��ư�� �����ֽñ� �ٶ��ϴ�.';
            DrawFormattedText(theWindow, double(msgtxt), 'center', 100, orange, [], [], [], 2);
            DrawFormattedText(theWindow, rating_types.prompts{ratetype}, 'center', 300, white, [], [], [], 2);

            [lb, rb, start_center] = draw_scale(exp_scale{1});

            [x,~,button] = GetMouse(theWindow);
            if x < lb; x = lb; elseif x > rb; x = rb; end
            if button(1); while button(1); [~,~,button] = GetMouse(theWindow); end; break; end
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end

            Screen('DrawLine', theWindow, orange, x, H/2, x, H/2+scale_W, 6);
            Screen('Flip', theWindow);
        end

        % Freeze the screen 0.5 second with red line
        Screen('DrawLine', theWindow, red, x, H/2, x, H/2+scale_W, 6);
        Screen('Flip', theWindow);

        freeze_t = GetSecs;
        while true
            freeze_cur_t = GetSecs;
            if freeze_cur_t - freeze_t > 0.5
                break
            end
        end


    case 'RATE2'

        msgtxt = sprintf(['�� �Կ��� %d ��° ������ ''���� ��'' �Դϴ�.\n\n', ...
            '�̹� ���������� �Կ��� �ռ�, ���� ������ ���� ��ü Ȱ���� �����ϰ� �˴ϴ�.\n', ...
            '�����ڴ� �����ڴ��� ��ü Ȱ���� �� ����Ǿ�����, ������ ����� ��ȭ�� �ִ��� Ȯ���մϴ�.'], runmarker);
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end

        msgtxt = ['��ü Ȱ���� ��ĳ�� �ۿ��� ����� ���, �Ӹ� ��ġ �����Կ��� �����մϴ�.'];
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end

        msgtxt = '+';
        DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end


    case 'RATE3'

        msgtxt = sprintf(['�� �Կ��� %d ��° ������ ''���� ��'' �Դϴ�.\n\n', ...
            '���� �������� ������ ���⸦ ���������� �����ֽø� �˴ϴ�.'], runmarker);
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end


    case 'PREP'

        msgtxt = sprintf(['�� �Կ��� %d ��° ������ ''���ϱ� �غ�'' �Դϴ�.\n\n', ...
            '�̹� ���������� �ٷ� ������ ���� ''���ϱ�'' �������� ���õ� �� ���� ������ ���ð� �˴ϴ�. \n', ...
            '�����ڴ��� �� ������ ���� �������� �亯���� ���� ������ �����ϸ鼭 ���� ������ �غ����ֽñ� �ٶ��ϴ�.'], runmarker);
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end
        
        if ~ismember(visitmarker, wh_tx)

            msgtxt = '������ ���� �������Դϴ�:';
            DrawFormattedText(theWindow, double(msgtxt), 'center', 100, orange, [], [], [], 2);
            msgtxt = ['1. ���� �����ϴ� �뷡�� � ���ΰ���? �� �뷡�� ������ ���� �Ǵ� ��ȭ�� �ִٸ� � ���ΰ���?\n\n', ...
                '2. �ֱٿ� �ô� ���� �� ���� ��ﳪ�� ���� �ֳ���? �� ������ ���� � ������ ������ ��̳���?'];
            DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
            Screen('Flip', theWindow);
            while true % Space
                [~,~,keyCode_E] = KbCheck(Exp_key);
                if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
                if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
            end
    
            msgtxt = ['��� �������� ������ ������ ���� ������, �亯�� ������ ������ ���������� ���谡 ��� �������ϴ�.\n', ...
                '�����ڴԲ��� ������ ���� �������� �Ǵ� ����Ǵ� � ���̵� �����ϴ�.\n\n', ...
                '�����ڴ��� �亯�Ͻ� ������ ���� ���� �������θ� ���Ǹ�, ö���� �͸�ȭ �� ��� ������ ����˴ϴ�.'];
            DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
            Screen('Flip', theWindow);
            while true % Space
                [~,~,keyCode_E] = KbCheck(Exp_key);
                if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
                if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
            end
    
            msgtxt = ['�� ���� �� �غ��ϼž� �ϴ� �亯�� ���̴� �� 5�о��Դϴ�.\n\n', ...
                '�亯�� �ʹ� ���� ������ ��쿡�� �߰� ������ �帱 �� �ֽ��ϴ�.'];
            DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
            Screen('Flip', theWindow);
            while true % Space
                [~,~,keyCode_E] = KbCheck(Exp_key);
                if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
                if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
            end

            msgtxt = ['���� �ʹ� ����ϱ� ����ϰų� �ڽſ��� ���� �ش���� �ʴ� �����̶��\n', ...
                '������ ��ư�� ������ ���� �������� ��ü�Ͻ� �� �ֽ��ϴ�.\n', ...
                '(���� ��ư�� ������ ���� �������� ���ư��ϴ�).\n\n', ...
                '�ٸ�, �غ�� ��ü ���� ������ �����Ǿ� �����Ƿ�, ���� ��ü �ÿ��� ������ ������ �ֽñ� �ٶ��ϴ�.'];
            DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
            Screen('Flip', theWindow);
            while true % Space
                [~,~,keyCode_E] = KbCheck(Exp_key);
                if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
                if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
            end

        elseif ismember(visitmarker, wh_tx)

            msgtxt = ['��� �������� ������ ������ ���� ������, �亯�� ������ ������ ���������� ���谡 ��� �������ϴ�.\n', ...
                '�����ڴԲ��� ������ ���� �������� �Ǵ� ����Ǵ� � ���̵� �����ϴ�.\n\n', ...
                '�����ڴ��� �亯�Ͻ� ������ ���� ���� �������θ� ���Ǹ�, ö���� �͸�ȭ �� ��� ������ ����˴ϴ�.'];
            DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
            Screen('Flip', theWindow);
            while true % Space
                [~,~,keyCode_E] = KbCheck(Exp_key);
                if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
                if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
            end

            msgtxt = ['������ ���õ� 2���� ������ 10, 20, 21-30ȸ�� �Կ����� ��� �����ϰ� �ݺ��� ���Դϴ�.\n', ...
                '����, ������ ����Ͻ� ���� ���� ������ ������ �������� ������ֽñ� �ٶ��ϴ�.\n\n', ...
                '�� ���� �� �غ��ϼž� �ϴ� �亯�� ���̴� �� 5�о��̸�, �ʹ� ���� ������ ��� �߰� ������ �帱 �� ������,\n', ...
                '�̹� ������ ���� ��� �ݺ��Ǵ� �߿��� �����̹Ƿ� �������̸� ����� ������ �亯�� �غ����ֽñ� ��Ź�帳�ϴ�.'];
            DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
            Screen('Flip', theWindow);
            while true % Space
                [~,~,keyCode_E] = KbCheck(Exp_key);
                if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
                if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
            end

        end
        

    case {'SPEAK1', 'SPEAK2'}

        msgtxt = sprintf(['�� �Կ��� %d ��° ������ ''���ϱ�'' �Դϴ�.\n\n', ...
            '�����ڴ��� ���� �������� ���õǾ��� %s ��° ������ ���� ����ϰ� �����Ӱ� �亯�Ͻø� �˴ϴ�.\n', ...
            '�亯 �ð��� �� 5���̸�, �亯�� �ʹ� ���� ������ ��쿡�� �����ڰ� �߰� ������ �帱 �� �ֽ��ϴ�.\n\n', ...
            '�����Ͻ� ������ ������ ������ �� ū ��Ҹ��� �������� ��Ȯ�ϰ� ���� ��Ź�帮��,\n', ...
            '�����Ͻô� ���� �Ӹ��� �����̽��� �ʵ��� ���� �ٶ��ϴ�.'], runmarker, S.type(end));
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end


    case 'LISTEN1'
        
        msgtxt = sprintf(['�� �Կ��� %d ��° ������ ''���'' �Դϴ�.\n\n', ...
            '�����ڴ��� ���� �������� %s ��° ������ ���� �亯�ϼ̴� ������ ������ ��� �˴ϴ�.\n', ...
            '�Կ� ����, ������ ���� ��Ÿ���� ȭ�� �߾��� + ǥ�ø� �����ϸ鼭 ������ �����ø� �˴ϴ�.'], runmarker, S.type(end));
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end

        msgtxt = '+';
        DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
        Screen('Flip', theWindow);
        temp_samp_file = which('funk.wav');
        [temp_samp_data, temp_samp_freq] = psychwavread(temp_samp_file);
        if size(temp_samp_data,2) == 1 && padev_output.NrOutputChannels == 2
            temp_samp_data = repmat(temp_samp_data, 1, 2);
        end
        pahandle = PsychPortAudio('Open', padev_output.DeviceIndex, 1, 0, temp_samp_freq, padev_output.NrOutputChannels);
        PsychPortAudio('FillBuffer', pahandle, temp_samp_data.');
        PsychPortAudio('Start', pahandle, 0, 0, 1);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end
        PsychPortAudio('Stop', pahandle);
        PsychPortAudio('Close');


    case 'LISTEN2'

        msgtxt = sprintf(['�� �Կ��� %d ��° ������ ''���'' �Դϴ�.\n\n', ...
            '�����ڴ��� ���� �������� %s ��° ������ ���� �亯�ϼ̴� ������ ������ ��� �˴ϴ�.\n', ...
            '�Կ� ����, ȭ�� �߾��� + ǥ�ø� �����ϸ鼭 ������ �����ø� �˴ϴ�.'], runmarker, S.type(end));
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end


    case 'RELISTEN1'

        msgtxt = sprintf(['�̹� ������ ''�ٽõ�� �� ��'' �Դϴ�.\n\n', ...
            '�����ڴ��� ���� ����ó��, %s ��° ������ ���� �亯 ���� ������ �ٽ� �� �� �� ��� �˴ϴ�.\n', ...
            '���� �������� ������ �����鼭 ������ ������ �ٽ� ���ø��� ô�� ���� ���������� ���Ͻø� �˴ϴ�.\n', ...
            '(���� ����� ������ �� �������� �ʴ´ٸ�, ���� �������� ������ �������� ���մϴ�)'], S.type(end));
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end

        % Explain scale with visualization
        msgtxt = '�ٽõ�� �� �� ����';
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, orange, [], [], [], 2);
        ratetype = strcmp(rating_types.alltypes, exp_scale{3});
        DrawFormattedText(theWindow, rating_types.prompts{ratetype}, 'center', 300, white, [], [], [], 2);
        [lb, rb, start_center] = draw_scale(exp_scale{3});
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end

        % Initial position
        if start_center
            SetMouse((rb+lb)/2,H/2); % set mouse at the center
        else
            SetMouse(lb,H/2); % set mouse at the left
        end

        temp_samp_file = which('funk.wav');
        [temp_samp_data, temp_samp_freq] = psychwavread(temp_samp_file);
        if size(temp_samp_data,2) == 1 && padev_output.NrOutputChannels == 2
            temp_samp_data = repmat(temp_samp_data, 1, 2);
        end
        pahandle = PsychPortAudio('Open', padev_output.DeviceIndex, 1, 0, temp_samp_freq, padev_output.NrOutputChannels);
        PsychPortAudio('FillBuffer', pahandle, temp_samp_data.');
        PsychPortAudio('Start', pahandle, 0, 0, 1);

        % Get ratings
        while true % Space
            msgtxt = '�����ڴ��� ����� �� ����� �����Ͻ� ��, ������ ������ ��ư�� �����ֽñ� �ٶ��ϴ�.';
            DrawFormattedText(theWindow, double(msgtxt), 'center', 100, orange, [], [], [], 2);
            DrawFormattedText(theWindow, rating_types.prompts{ratetype}, 'center', 300, white, [], [], [], 2);

            [lb, rb, start_center] = draw_scale(exp_scale{3});

            [x,~,button] = GetMouse(theWindow);
            if x < lb; x = lb; elseif x > rb; x = rb; end
            if button(1); while button(1); [~,~,button] = GetMouse(theWindow); end; break; end
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end

            Screen('DrawLine', theWindow, orange, x, H/2, x, H/2+scale_W, 6);
            Screen('Flip', theWindow);
        end

        % Freeze the screen 0.5 second with red line
        Screen('DrawLine', theWindow, red, x, H/2, x, H/2+scale_W, 6);
        Screen('Flip', theWindow);

        freeze_t = GetSecs;
        while true
            freeze_cur_t = GetSecs;
            if freeze_cur_t - freeze_t > 0.5
                break
            end
        end

        PsychPortAudio('Stop', pahandle);
        PsychPortAudio('Close');


    case 'RELISTEN2'

        msgtxt = sprintf(['�̹� ������ ''�ٽõ�� �� ��'' �Դϴ�.\n\n', ...
            '�����ڴ��� ���� ����ó��, %s ��° ������ ���� �亯 ���� ������ �ٽ� �� �� �� ��� �˴ϴ�.\n', ...
            '���� �������� ������ �����鼭 �������� ������ �ٽ� ���ø��� ô�� ���� �����ֽø� �˴ϴ�.'], S.type(end));
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end


end


%% SETUP: Sound Setting (for SPEAK, LISTEN, and RELISTEN)

if ismember(S.type, {'SPEAK1', 'SPEAK2'})

    pahandle = PsychPortAudio('Open', padev_input.DeviceIndex, 2, 1, [], padev_input.NrInputChannels);
    data.dat.recording_struct = PsychPortAudio('GetStatus', pahandle);
    PsychPortAudio('GetAudioData', pahandle, S.dur + 10); % actual = S.dur + dummy 6 secs
    PsychPortAudio('Start', pahandle, 0, 0, 1);
    PsychPortAudio('Stop', pahandle);
    PsychPortAudio('GetAudioData', pahandle);

elseif ismember(S.type, {'LISTEN1', 'LISTEN2', 'RELISTEN1', 'RELISTEN2'})

    prevdatafile = fullfile(savedir, sprintf('%s_%s_sub_%.3d_visit_%.2d_run_%.2d.mat', subjtime, subjID, subjnum, visitmarker, find(strcmp(runs, ['SPEAK' S.type(end)]))));
    prevdata = load(prevdatafile, 'data');
    [~, prevdata_name] = fileparts(strrep(strrep(prevdata.data.dat.recording_file,'\',filesep),'/',filesep));
    data.dat.recording_proc_file = fullfile(basedir, 'Recording', [prevdata_name '_p.wav']);
    [rec_proc_data, rec_proc_freq] = psychwavread(data.dat.recording_proc_file);
    if size(rec_proc_data,2) == 1 && padev_output.NrOutputChannels == 2
        rec_proc_data = repmat(rec_proc_data, 1, 2);
    end
    rec_proc_data = rec_proc_data(ceil(rec_proc_freq * prevdata.data.dat.recording_dummydur)+1 : end, :); % Remove dummy 6 seconds
    pahandle = PsychPortAudio('Open', padev_output.DeviceIndex, 1, 0, rec_proc_freq, padev_output.NrOutputChannels);
    PsychPortAudio('FillBuffer', pahandle, rec_proc_data.');
    data.dat.playing_preptime = PsychPortAudio('Start', pahandle, 1, 0, 1, GetSecs);
    PsychPortAudio('Stop', pahandle);

end


%% Main : Ready for scan

if ~ismember(S.type, {'RELISTEN1', 'RELISTEN2'})

    msgtxt = ['���ݺ��� �� ������ ���۵˴ϴ�.\n\n', ...
        '����: �Կ� �� �Ӹ��� �����̰ų� �ῡ ���� �ʵ��� �������ֽñ� �ٶ��ϴ�!!!'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, orange, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end
    
    msgtxt = ['�����ڴ� ��� ���� �� �����ڴ��� �غ� �Ϸ�Ǿ����� Ȯ���ϱ� �ٶ��ϴ�.\n', ...
        '�غ� �Ϸ�Ǹ� �����ڴ� SPACE Ű�� �����ֽñ� �ٶ��ϴ�.'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    msgtxt = ['�ڱ��� ���� �Կ��� �ʿ��� ��� �����մϴ�.'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    msgtxt = '+';
    DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

elseif ismember(S.type, {'RELISTEN1', 'RELISTEN2'})

    msgtxt = ['���ݺ��� �� ������ ���۵˴ϴ�.\n\n', ...
        '�����ڴԲ����� ��� �غ� �Ϸ��Ͻ� �� SPACE Ű�� �����ֽñ� �ٶ��ϴ�.'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

end

if ~strcmp(S.type, 'PREP')
    %% MAIN : Sync (S key)
    
    if ~ismember(S.type, {'RELISTEN1', 'RELISTEN2'})
        msgtxt = '��ĵ�� �����մϴ�. (S Ű)';
        DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true
            [~,~,keyCode_S] = KbCheck(Scan_key);
            if keyCode_S(KbName('s')); break; end
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end
    end
    
    
    %% MAIN : Disdaq (4 + 4 + 2 = 10 secs)
    
    % 4 secs : scanning...
    start_t = GetSecs;
    data.dat.runscan_starttime = start_t;

    if ismember(S.type, {'SPEAK1', 'SPEAK2'})
        data.dat.recording_preptime = PsychPortAudio('Start', pahandle, 0, 0, 1);
    end
    
    msgtxt = '�����ϴ� ��...';
    DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
    Screen('Flip', theWindow);

    while true
        cur_t = GetSecs;
        if cur_t - start_t >= 4
            break
        end
    end
    
    if ismember(S.type, {'SPEAK1', 'SPEAK2'})
        PsychPortAudio('GetAudioData', pahandle);
        data.dat.recording_starttime = GetSecs;
    end
    
    % 4 secs : blank
    Screen(theWindow, 'FillRect', bgcolor, window_rect);
    Screen('Flip', theWindow);
    while true
        cur_t = GetSecs;
        if cur_t - start_t >= 8
            break
        end
    end
    
    % 2 secs : BIOPAC
    if USE_BIOPAC && ~ismember(S.type, {'PREP', 'RELISTEN1', 'RELISTEN2'})
        BIOPAC_trigger(ljHandle, biopac_channel, 'on');
    end
    while true
        cur_t = GetSecs;
        if cur_t - start_t >= 10
            break
        end
    end
    if USE_BIOPAC && ~ismember(S.type, {'PREP', 'RELISTEN1', 'RELISTEN2'})
        BIOPAC_trigger(ljHandle, biopac_channel, 'off');
    end

elseif strcmp(S.type, 'PREP')

    msgtxt = '��ĵ�� �����մϴ�. (Space Ű)';
    DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

end


%% MAIN : Experiment

start_t = GetSecs;
data.dat.experiment_starttime = start_t;
if ismember(S.type, {'SPEAK1', 'SPEAK2'})
    data.dat.recording_dummydur = data.dat.experiment_starttime - data.dat.recording_starttime;
    fprintf('\n\n*************\n\nDUMMY DURATION: %.3f\n\n*************\n\n', data.dat.recording_dummydur);
end

switch S.type


    case {'REST'}

        msgtxt = '+';
        DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            cur_t = GetSecs;
            if cur_t - start_t >= S.dur
                break
            end
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end


    case {'RATE1', 'RATE2', 'RATE3'}

        % Basic setting
        rec_i = 0;
        ratetype = strcmp(rating_types.alltypes, main_scale{1});

        [lb, rb, start_center] = draw_scale(main_scale{1}); % Getting information
        Screen('FillRect', theWindow, bgcolor, window_rect); % clear the screen

        % Initial position
        if start_center
            SetMouse((rb+lb)/2,H/2); % set mouse at the center
        else
            SetMouse(lb,H/2); % set mouse at the left
        end


        % Get ratings
        time_fromstart = NaN(40000,1); % ~11 min given 60Hz flip freq
        cont_rating = NaN(40000,1); % ~11 min given 60Hz flip freq
        data.dat.rating_starttime = GetSecs;
        while true
            DrawFormattedText(theWindow, rating_types.prompts{ratetype}, 'center', 200, orange, [], [], [], 2);

            [lb, rb, start_center] = draw_scale(main_scale{1});

            [x,~,~] = GetMouse(theWindow);
            if x < lb; x = lb; elseif x > rb; x = rb; end
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end

            rec_i = rec_i + 1;
            cur_t = GetSecs;
            time_fromstart(rec_i) = cur_t-start_t;
            cont_rating(rec_i) = (x-lb)./(rb-lb);

            if cur_t - start_t >= S.dur
                break
            end

            Screen('DrawLine', theWindow, orange, x, H/2, x, H/2+scale_W, 6);
            Screen('Flip', theWindow);
        end
        data.dat.time_fromstart = time_fromstart(1:rec_i);
        data.dat.cont_rating = cont_rating(1:rec_i);


    case {'PREP'}

        if ~ismember(visitmarker, wh_tx)
            data.dat.questiontype = 'semantic';
            data.dat.questionmarker_before = questionmarker;
            while true
                change_q = false;
                msgtxt = sprintf('1. %s\n\n2. %s', question_all.semantic{questionmarker,:});
                DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
                Screen('Flip', theWindow);
                while true % Space
                    [~,~,button] = GetMouse(theWindow);
                    if button(2) || button(3) % Rt click then next question
                        while button(2) || button(3); [~,~,button] = GetMouse(theWindow); end
                        questionmarker = min(questionmarker + 1, visitmarker + 5);
                        change_q = true;
                        break
                    elseif button(1) % Lt click then previous question
                        while button(1); [~,~,button] = GetMouse(theWindow); end
                        questionmarker = max(questionmarker - 1, data.dat.questionmarker_before);
                        change_q = true;
                        break
                    end
                    cur_t = GetSecs;
                    if cur_t - start_t >= S.dur
                        break
                    end
                    [~,~,keyCode_E] = KbCheck(Exp_key);
                    if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
                end
                if ~change_q; break; end
            end
            data.dat.questionmarker_after = questionmarker;
            data.dat.question = msgtxt;
        elseif ismember(visitmarker, wh_tx)
            data.dat.questiontype = 'therapeutic';
            msgtxt = sprintf('1. %s\n\n2. %s', question_all.therapeutic{:});
            DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
            Screen('Flip', theWindow);
            while true % Space
                cur_t = GetSecs;
                if cur_t - start_t >= S.dur
                    break
                end
            end
            data.dat.question = msgtxt;
        end

        msgtxt = ['���ϱ� �غ� �������ϴ�.\n\n', ...
            '��ĵ�� �Ϸ�Ǹ� �����ڴ� SPACE Ű�� �����ֽñ� �ٶ��ϴ�.'];
        DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end


    case {'SPEAK1', 'SPEAK2'}
        
        if ~ismember(visitmarker, wh_tx)
            data.dat.questiontype = 'semantic';
            data.dat.questionmarker = questionmarker;
            msgtxt = question_all.semantic{questionmarker, str2double(S.type(end))};
        elseif ismember(visitmarker, wh_tx)
            data.dat.questiontype = 'therapeutic';
            msgtxt = question_all.therapeutic{str2double(S.type(end))};
        end

        DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
        Screen('Flip', theWindow);
        % Recording
        data.dat.recording_speakstarttime = GetSecs;
        while true
            cur_t = GetSecs;
            if cur_t - start_t >= S.dur
                break
            end
        end
        rec_data = PsychPortAudio('GetAudioData', pahandle);
        rec_data = rec_data(1,:); % Only Lt channel - less noisy
        data.dat.recording_endtime = GetSecs;
        data.dat.question = msgtxt;
        PsychPortAudio('Stop', pahandle);
        PsychPortAudio('Close');
        data.dat.recording_file = data.datafile;
        data.dat.recording_file = strrep(data.dat.recording_file, savedir, fullfile(basedir, 'Recording'));
        data.dat.recording_file = strrep(data.dat.recording_file, '.mat', '_record.wav');
        psychwavwrite(rec_data.', data.dat.recording_struct.SampleRate, data.dat.recording_file);


    case {'LISTEN1', 'LISTEN2'}
  
        msgtxt = '+';
        DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
        Screen('Flip', theWindow);
        % Playing
        data.dat.playing_starttime = PsychPortAudio('Start', pahandle, 1, 0, 1);
        while true % Space
            cur_t = GetSecs;
            if cur_t - start_t >= S.dur
                break
            end
        end
        PsychPortAudio('Stop', pahandle);
        data.dat.playing_endtime = GetSecs;
        PsychPortAudio('Close'); % This takes some time. Exact dur of recording is data.dat.time_fromstart(end).

    case {'RELISTEN1', 'RELISTEN2'}

        % Basic setting
        rec_i = 0;
        ratetype = strcmp(rating_types.alltypes, main_scale{2});

        [lb, rb, start_center] = draw_scale(main_scale{2}); % Getting information
        Screen('FillRect', theWindow, bgcolor, window_rect); % clear the screen

        % Initial position
        if start_center
            SetMouse((rb+lb)/2,H/2); % set mouse at the center
        else
            SetMouse(lb,H/2); % set mouse at the left
        end

        data.dat.playing_starttime = PsychPortAudio('Start', pahandle, 1, 0, 1);

        % Get ratings
        time_fromstart = NaN(20000,1); % ~5.5 min given 60Hz flip freq
        cont_rating = NaN(20000,1); % ~5.5 min given 60Hz flip freq
        data.dat.rating_starttime = GetSecs;
        while true
            DrawFormattedText(theWindow, rating_types.prompts{ratetype}, 'center', 200, orange, [], [], [], 2);

            [lb, rb, start_center] = draw_scale(main_scale{2});

            [x,~,~] = GetMouse(theWindow);
            if x < lb; x = lb; elseif x > rb; x = rb; end
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end

            rec_i = rec_i + 1;
            cur_t = GetSecs;
            time_fromstart(rec_i) = cur_t-start_t;
            cont_rating(rec_i) = (x-lb)./(rb-lb);

            if cur_t - start_t >= S.dur
                break
            end

            Screen('DrawLine', theWindow, orange, x, H/2, x, H/2+scale_W, 6);
            Screen('Flip', theWindow);
        end
        data.dat.time_fromstart = time_fromstart(1:rec_i);
        data.dat.cont_rating = cont_rating(1:rec_i);

        PsychPortAudio('Stop', pahandle);
        data.dat.playing_endtime = GetSecs;
        PsychPortAudio('Close');

end

cur_t = GetSecs;
data.dat.experiment_endtime = cur_t;
data.dat.experiment_total_dur = cur_t - start_t;

if USE_BIOPAC && ~ismember(S.type, {'PREP', 'RELISTEN1', 'RELISTEN2'})
    BIOPAC_trigger(ljHandle, biopac_channel, 'on');
    BIOPAC_trigger(ljHandle, biopac_channel, 'off');
end

save(data.datafile, '-append', 'data');


%% MAIN : Postrun questionnaire

all_start_t = GetSecs;
data.dat.postrun_rating_timestamp = all_start_t;
ratestim = strcmp(rating_types.postallstims, S.type);
scales = rating_types.postalltypes{ratestim};

Screen(theWindow, 'FillRect', bgcolor, window_rect);
Screen('Flip', theWindow); % clear screen

% Going through each scale
for scale_i = 1:numel(scales)
    
    % First introduction
    if scale_i == 1
        
        msgtxt = [num2str(runmarker) '��° ������ �������ϴ�.\n��� �� �������� ���õ� ���Դϴ�. �����ںв����� ��ٷ��ֽñ� �ٶ��ϴ�.'];
        DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
        Screen('Flip', theWindow);
        
        start_t = GetSecs;
        while true
            cur_t = GetSecs;
            if cur_t - start_t >= postrun_start_t
                break
            end
        end
        
    end
    
    % Parse scales and basic setting
    scale = scales{scale_i};
    
    [lb, rb, start_center] = draw_scale(scale);
    Screen(theWindow, 'FillRect', bgcolor, window_rect);
    
    start_t = GetSecs;
    data.dat = setfield(data.dat, sprintf('%s_timestamp', scale), start_t);
    
    rec_i = 0;
    ratetype = strcmp(rating_types.alltypes, scale);
    
    % Initial position
    if start_center
        SetMouse((rb+lb)/2,H/2); % set mouse at the center
    else
        SetMouse(lb,H/2); % set mouse at the left
    end
    
    % Get ratings
    while true % Button
        DrawFormattedText(theWindow, rating_types.prompts{ratetype}, 'center', 200, white, [], [], [], 2);

        [lb, rb, start_center] = draw_scale(scale);
        
        [x,~,button] = GetMouse(theWindow);
        if x < lb; x = lb; elseif x > rb; x = rb; end
        if button(1); while button(1); [~,~,button] = GetMouse(theWindow); end; break; end
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        
        Screen('DrawLine', theWindow, orange, x, H/2, x, H/2+scale_W, 6);
        Screen('Flip', theWindow);
    end
    
    end_t = GetSecs;
    data.dat = setfield(data.dat, sprintf('%s_rating', scale), (x-lb)./(rb-lb));
    data.dat = setfield(data.dat, sprintf('%s_RT', scale), end_t - start_t);

    % Freeze the screen 0.5 second with red line
    Screen('DrawLine', theWindow, red, x, H/2, x, H/2+scale_W, 6);
    Screen('Flip', theWindow);

    freeze_t = GetSecs;
    while true
        freeze_cur_t = GetSecs;
        if freeze_cur_t - freeze_t > 0.5
            break
        end
    end
    
    if scale_i == numel(scales)

        msgtxt = '������ �������ϴ�.';
        DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
        Screen('Flip', theWindow);
        
        start_t = GetSecs;
        while true
            cur_t = GetSecs;
            if cur_t - start_t >= postrun_end_t
                break
            end
        end
        
    end
    
end

all_end_t = GetSecs;
data.dat.postrun_total_RT = all_end_t - all_start_t;

save(data.datafile, '-append', 'data');


%% Closing screen

msgtxt = [num2str(runmarker) '��° ������ �������ϴ�.\n������ ��ġ����, �����ڴ� SPACE Ű�� �����ֽñ� �ٶ��ϴ�.'];
DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
Screen('Flip', theWindow);
while true % Space
    [~,~,keyCode_E] = KbCheck(Exp_key);
    if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
    if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
end
Screen('CloseAll');


%% Update markers and finish experiment

if runmarker < size(marker_mat,3)
    runmarker = runmarker + 1;
elseif runmarker == size(marker_mat,3)
    visitmarker = visitmarker + 1;
    runmarker = 1;
end

marker_mat(subjnum, :, :) = false;
if visitmarker <= size(marker_mat,2)
    marker_mat(subjnum, visitmarker, runmarker) = true;
end

if strcmp(S.type, 'SPEAK2') && ~ismember(visitmarker, wh_tx)
    questionmarker = questionmarker + 1;
end
qmarker_mat(subjnum, :) = false;
qmarker_mat(subjnum, questionmarker) = true;

save(markerfile, '-append', 'marker_mat', 'qmarker_mat');

disp('Done');
