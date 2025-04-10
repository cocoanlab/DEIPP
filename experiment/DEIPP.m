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
        inputdev.DeviceName = 'JaeJoong Lee의 AirPods'; %'내장 마이크'
        outputdev.HostAudioAPIName = 'Core Audio';
        outputdev.DeviceName = 'JaeJoong Lee의 AirPods'; %'내장 출력'
    case 'Cocoanui-iMac-2.local'
        basedir = '/Users/deipp/Dropbox/Paperwork/DEIPP_sync/experiment/DEIPP_ex';
        inputdev.HostAudioAPIName = 'Core Audio';
        inputdev.DeviceName = '내장 마이크';
        outputdev.HostAudioAPIName = 'Core Audio';
        outputdev.DeviceName = '내장 출력';
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
[~, ~, wordrect1, ~] = DrawFormattedText(theWindow, double('코'), lb1-30, H/2+scale_W+40, bgcolor);
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

    msgtxt = ['안녕하세요, 참가자님. ''뇌 기능 영상을 통한 만성 통증 연구''에 참가해주셔서 감사합니다.\n\n', ...
        '본 촬영은 다음과 같은 총 4개의 과제로 구성되어 있으며, 약 1시간 30분이 소요됩니다.'];
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

    msgtxt = ['참가자님은 촬영에 앞서 최대한 편한 자세를 취하신 후, 실험이 끝날 때까지 유지해주시기 바라며,\n' ...
        '머리를 움직이거나 잠에 들지 않도록 각별히 유의해주시기 바랍니다.'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    msgtxt = ['먼저, 참가자님의 머리 위치 파악을 위한 예비 촬영을 진행하겠습니다.\n\n', ...
        '참가자님은 다음과 같이 나타나는 화면 중앙의 + 표시를 응시하면서 편안히 계시면 됩니다.'];
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

    msgtxt = ['다음으로는, 촬영 중 소음 감소 테스트입니다.\n\n', ...
        '참가자님은 화면 중앙의 + 표시를 응시하면서, 실험자의 안내를 기다려주시기 바랍니다.'];
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
    msgtxt = '스캐너의 소음이 충분히 줄어들었는지, 실험자와 참가자의 소리 크기가 적절한 지 확인합니다.';
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    % Finish testing
    msgtxt = ['테스트가 완료되었습니다. 실험자는 테스트 스캔을 종료합니다.\n\n', ...
        '지금부터 본 실험이 시작됩니다.'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

elseif strcmp(S.type, 'SPEAK1') && check_recording

    msgtxt = ['촬영에 앞서, 소리 입력-출력 테스트를 진행하겠습니다.\n\n', ...
        '참가자님은 화면 중앙의 + 표시를 응시하면서, 실험자의 안내를 기다려주시기 바랍니다.'];
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
    msgtxt = ['실험자와 참가자의 사이의 대화 녹음을 테스트합니다.\n\n', ...
        '화면에 + 표시가 나타나면 실험자의 질문에 대답해주세요!']; % ask name and then age
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
    msgtxt = ['녹음 파일 저장 및 스캐너 소음 제거 중입니다...\n\n', ...
        '참가자님은 조금만 기다려 주시기 바랍니다.'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    % Check sound volume of processed recording data (if needed, adjust the Laptop or Line 1 volume)
    msgtxt = ['완료되었습니다. 녹음 파일을 재생합니다.\n\n', ...
        '참가자님은 녹음을 듣고 소리 크기가 적절한 지 확인합니다.'];
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
    msgtxt = ['테스트가 완료되었습니다. 실험자는 테스트 스캔을 종료합니다.\n\n', ...
        '지금부터 본 실험이 시작됩니다.'];
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

        msgtxt = sprintf(['본 촬영의 %d 번째 순서는 ''휴지기 촬영'' 입니다.\n\n', ...
            '휴지기 촬영 동안, 참가자님은 다음과 같이 나타나는 화면 중앙의 + 표시를 응시하면서 편안히 계시면 됩니다.'], runmarker);
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
        msgtxt = '각 촬영 세션이 끝날 때마다, 간단한 질문 및 평가가 이루어집니다.';
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
            msgtxt = '참가자님은 충분히 평가 방법을 연습하신 후, 연습이 끝나면 버튼을 눌러주시기 바랍니다.';
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

        msgtxt = sprintf(['본 촬영의 %d 번째 순서는 ''통증 평가'' 입니다.\n\n', ...
            '촬영 동안 참가자님은 다음과 같은 척도 위에 현재 느껴지는 통증의 세기를 지속적으로 평가해주시면 됩니다.'], runmarker);
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end

        % Explain scale with visualization
        msgtxt = '통증 평가 과제';
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
            msgtxt = '참가자님은 충분히 평가 방법을 연습하신 후, 연습이 끝나면 버튼을 눌러주시기 바랍니다.';
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

        msgtxt = sprintf(['본 촬영의 %d 번째 순서는 ''통증 평가'' 입니다.\n\n', ...
            '이번 순서에서는 촬영에 앞서, 먼저 사전에 정한 신체 활동을 시행하게 됩니다.\n', ...
            '실험자는 참가자님의 신체 활동이 잘 진행되었는지, 통증에 충분한 변화가 있는지 확인합니다.'], runmarker);
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end

        msgtxt = ['신체 활동이 스캐너 밖에서 진행된 경우, 머리 위치 예비촬영을 진행합니다.'];
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

        msgtxt = sprintf(['본 촬영의 %d 번째 순서는 ''통증 평가'' 입니다.\n\n', ...
            '현재 느껴지는 통증의 세기를 지속적으로 평가해주시면 됩니다.'], runmarker);
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end


    case 'PREP'

        msgtxt = sprintf(['본 촬영의 %d 번째 순서는 ''말하기 준비'' 입니다.\n\n', ...
            '이번 순서에서는 바로 다음에 있을 ''말하기'' 순서에서 제시될 두 가지 질문을 보시게 됩니다. \n', ...
            '참가자님은 각 질문에 대해 떠오르는 답변들을 마음 속으로 생각하면서 다음 순서를 준비해주시기 바랍니다.'], runmarker);
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end
        
        if ~ismember(visitmarker, wh_tx)

            msgtxt = '다음은 예시 질문들입니다:';
            DrawFormattedText(theWindow, double(msgtxt), 'center', 100, orange, [], [], [], 2);
            msgtxt = ['1. 가장 좋아하는 노래가 어떤 것인가요? 그 노래와 연관된 경험 또는 일화가 있다면 어떤 것인가요?\n\n', ...
                '2. 최근에 봤던 뉴스 중 가장 기억나는 것이 있나요? 그 뉴스를 보고 어떤 생각과 감정이 드셨나요?'];
            DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
            Screen('Flip', theWindow);
            while true % Space
                [~,~,keyCode_E] = KbCheck(Exp_key);
                if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
                if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
            end
    
            msgtxt = ['모든 질문에는 정답이 정해져 있지 않으며, 답변의 내용이 질문과 직접적으로 관계가 없어도 괜찮습니다.\n', ...
                '참가자님께서 질문을 보고 생각나는 또는 연상되는 어떤 것이든 좋습니다.\n\n', ...
                '참가자님이 답변하신 내용은 오직 연구 목적으로만 사용되며, 철저한 익명화 및 비밀 유지가 보장됩니다.'];
            DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
            Screen('Flip', theWindow);
            while true % Space
                [~,~,keyCode_E] = KbCheck(Exp_key);
                if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
                if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
            end
    
            msgtxt = ['각 질문 당 준비하셔야 하는 답변의 길이는 총 5분씩입니다.\n\n', ...
                '답변이 너무 일찍 끝나는 경우에는 추가 질문을 드릴 수 있습니다.'];
            DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
            Screen('Flip', theWindow);
            while true % Space
                [~,~,keyCode_E] = KbCheck(Exp_key);
                if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
                if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
            end

            msgtxt = ['만약 너무 대답하기 곤란하거나 자신에게 전혀 해당되지 않는 질문이라면\n', ...
                '오른쪽 버튼을 눌러서 다음 질문으로 교체하실 수 있습니다.\n', ...
                '(왼쪽 버튼을 누르면 이전 질문으로 돌아갑니다).\n\n', ...
                '다만, 준비된 전체 질문 갯수는 한정되어 있으므로, 질문 교체 시에는 신중히 결정해 주시길 바랍니다.'];
            DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
            Screen('Flip', theWindow);
            while true % Space
                [~,~,keyCode_E] = KbCheck(Exp_key);
                if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
                if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
            end

        elseif ismember(visitmarker, wh_tx)

            msgtxt = ['모든 질문에는 정답이 정해져 있지 않으며, 답변의 내용이 질문과 직접적으로 관계가 없어도 괜찮습니다.\n', ...
                '참가자님께서 질문을 보고 생각나는 또는 연상되는 어떤 것이든 좋습니다.\n\n', ...
                '참가자님이 답변하신 내용은 오직 연구 목적으로만 사용되며, 철저한 익명화 및 비밀 유지가 보장됩니다.'];
            DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
            Screen('Flip', theWindow);
            while true % Space
                [~,~,keyCode_E] = KbCheck(Exp_key);
                if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
                if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
            end

            msgtxt = ['앞으로 제시될 2개의 질문은 10, 20, 21-30회차 촬영에서 계속 동일하게 반복될 것입니다.\n', ...
                '따라서, 질문에 대답하실 때는 현재 본인의 생각을 기준으로 대답해주시기 바랍니다.\n\n', ...
                '각 질문 당 준비하셔야 하는 답변의 길이는 총 5분씩이며, 너무 일찍 끝나는 경우 추가 질문을 드릴 수 있으나,\n', ...
                '이번 질문은 이후 계속 반복되는 중요한 질문이므로 가급적이면 충분한 길이의 답변을 준비해주시기 부탁드립니다.'];
            DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
            Screen('Flip', theWindow);
            while true % Space
                [~,~,keyCode_E] = KbCheck(Exp_key);
                if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
                if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
            end

        end
        

    case {'SPEAK1', 'SPEAK2'}

        msgtxt = sprintf(['본 촬영의 %d 번째 순서는 ''말하기'' 입니다.\n\n', ...
            '참가자님은 이전 순서에서 제시되었던 %s 번째 질문에 대해 편안하고 자유롭게 답변하시면 됩니다.\n', ...
            '답변 시간은 총 5분이며, 답변이 너무 일찍 끝나는 경우에는 실험자가 추가 질문을 드릴 수 있습니다.\n\n', ...
            '말씀하실 때에는 언제나 가능한 한 큰 목소리로 느리더라도 정확하게 말씀 부탁드리며,\n', ...
            '말씀하시는 동안 머리는 움직이시지 않도록 주의 바랍니다.'], runmarker, S.type(end));
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end


    case 'LISTEN1'
        
        msgtxt = sprintf(['본 촬영의 %d 번째 순서는 ''듣기'' 입니다.\n\n', ...
            '참가자님은 이전 순서에서 %s 번째 질문에 대해 답변하셨던 내용의 녹음을 듣게 됩니다.\n', ...
            '촬영 동안, 다음과 같이 나타나는 화면 중앙의 + 표시를 응시하면서 녹음을 들으시면 됩니다.'], runmarker, S.type(end));
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

        msgtxt = sprintf(['본 촬영의 %d 번째 순서는 ''듣기'' 입니다.\n\n', ...
            '참가자님은 이전 순서에서 %s 번째 질문에 대해 답변하셨던 내용의 녹음을 듣게 됩니다.\n', ...
            '촬영 동안, 화면 중앙의 + 표시를 응시하면서 녹음을 들으시면 됩니다.'], runmarker, S.type(end));
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end


    case 'RELISTEN1'

        msgtxt = sprintf(['이번 과제는 ''다시듣기 및 평가'' 입니다.\n\n', ...
            '참가자님은 이전 과제처럼, %s 번째 질문에 대한 답변 녹음 내용을 다시 한 번 더 듣게 됩니다.\n', ...
            '이전 과제에서 녹음을 들으면서 느꼈던 감정을 다시 떠올리며 척도 위에 지속적으로 평가하시면 됩니다.\n', ...
            '(만약 당시의 감정이 잘 떠오르지 않는다면, 현재 느껴지는 감정을 기준으로 평가합니다)'], S.type(end));
        DrawFormattedText(theWindow, double(msgtxt), 'center', 100, white, [], [], [], 2);
        Screen('Flip', theWindow);
        while true % Space
            [~,~,keyCode_E] = KbCheck(Exp_key);
            if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
            if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
        end

        % Explain scale with visualization
        msgtxt = '다시듣기 및 평가 과제';
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
            msgtxt = '참가자님은 충분히 평가 방법을 연습하신 후, 연습이 끝나면 버튼을 눌러주시기 바랍니다.';
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

        msgtxt = sprintf(['이번 과제는 ''다시듣기 및 평가'' 입니다.\n\n', ...
            '참가자님은 이전 과제처럼, %s 번째 질문에 대한 답변 녹음 내용을 다시 한 번 더 듣게 됩니다.\n', ...
            '이전 과제에서 녹음을 들으면서 느꼈었던 감정을 다시 떠올리며 척도 위에 평가해주시면 됩니다.'], S.type(end));
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

    msgtxt = ['지금부터 본 실험이 시작됩니다.\n\n', ...
        '주의: 촬영 중 머리를 움직이거나 잠에 들지 않도록 유의해주시기 바랍니다!!!'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 100, orange, [], [], [], 2);
    Screen('Flip', theWindow);
    while true % Space
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end
    
    msgtxt = ['실험자는 모든 세팅 및 참가자님의 준비가 완료되었는지 확인하기 바랍니다.\n', ...
        '준비가 완료되면 실험자는 SPACE 키를 눌러주시기 바랍니다.'];
    DrawFormattedText(theWindow, double(msgtxt), 'center', 'center', white, [], [], [], 2);
    Screen('Flip', theWindow);
    while true
        [~,~,keyCode_E] = KbCheck(Exp_key);
        if keyCode_E(KbName('space')); while keyCode_E(KbName('space')); [~,~,keyCode_E] = KbCheck(Exp_key); end; break; end
        if keyCode_E(KbName('q')); abort_experiment('manual'); break; end
    end

    msgtxt = ['자기장 보정 촬영이 필요한 경우 진행합니다.'];
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

    msgtxt = ['지금부터 본 실험이 시작됩니다.\n\n', ...
        '참가자님께서는 모든 준비를 완료하신 후 SPACE 키를 눌러주시기 바랍니다.'];
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
        msgtxt = '스캔을 시작합니다. (S 키)';
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
    
    msgtxt = '시작하는 중...';
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

    msgtxt = '스캔을 시작합니다. (Space 키)';
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

        msgtxt = ['말하기 준비가 끝났습니다.\n\n', ...
            '스캔이 완료되면 실험자는 SPACE 키를 눌러주시기 바랍니다.'];
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
        
        msgtxt = [num2str(runmarker) '번째 세션이 끝났습니다.\n잠시 후 질문들이 제시될 것입니다. 참가자분께서는 기다려주시기 바랍니다.'];
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

        msgtxt = '질문이 끝났습니다.';
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

msgtxt = [num2str(runmarker) '번째 세션이 끝났습니다.\n세션을 마치려면, 실험자는 SPACE 키를 눌러주시기 바랍니다.'];
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
