% provides the instructions and test questions for the Ultimatum Game and
% for the estimation of the probability task as well as the "exchange task"
% (i.e. the dictator game)
% Author: Michael Giffin, January, 2018
%------ inputs -------
clearvars;
% stuff the script will ask me for before running the experiment
sub_num = input('What is the subject number?\n');

age     = input('What is the subject''s age?\n');
gender     = input('What is the subject''s gender?\n');

% sub_num = 1;
% for now I'm just running it as condition 1
% conditions are here for counterbalancing
%
%           condition 1: social nonsocial social nonsocial; social filled nonsocial
%           framed
%
%           condition 2: nonsocial social nonsocial social; social filled
%           nonsocial framed
%
%           condition 3: social nonsocial social nonsocial; social framed nonsocial
%           filled
%
%           condition 4: nonsocial social nonsocial social; social framed
%           nonsocial filled
%
% condition = input('Condition (1:4): \n');
% condition = 1;
condition = mod(sub_num, 4); % condition 4 will be recorded as a 0 with mod function



% opponent images
if condition == 1
    socmat = [1 0 1 0];
    socpo  = 'filled';
    nonsocpo= 'frame';
elseif condition == 2
    socmat = [0 1 0 1];
    socpo  = 'filled';
    nonsocpo= 'frame';
elseif condition == 3
    socmat = [1 0 1 0];
    socpo  = 'frame';
    nonsocpo= 'filled';
elseif condition == 0
    socmat = [0 1 0 1];
    socpo  = 'frame';
    nonsocpo= 'filled';
end


% Clear the workspace and the screen
sca;
close all;

script_start = GetSecs;


%=========================================================================
% VERY IMPORTANT - MUST BE COMMENTED OUT FOR REAL EXPERIMENT
%
% this I am including because apparently there's an issue with proppper
% timing and Mac High Sierra, luckily I won't be using this computer for
% the data collection (in fact the computer that's being used in the
% scanner is a windows computer). Therefore, when writing and testing the
% code on my computer this can stay (unless I can get things to work when
% my computer is attached to an external monitor), but must must MUST be
% commented out in the real experiment to ensure sensitive timing

Screen('Preference', 'SkipSyncTests', 1);

%=========================================================================


% have function that shuffles up the block image directories on the first
% block so that all the images are randomized

% this is supposed to make the keyboard specific to the operating system
% I'm on and should generalize accross operating systems
KbName('UnifyKeyNames');


% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);

% Define black and white

% Define black, white and grey
white = WhiteIndex(screenNumber);
grey = white / 2;
black = BlackIndex(screenNumber);


% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);

center        = round([windowRect(3) windowRect(4)]/2);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');


% Query the maximum priority level
topPriorityLevel = MaxPriority(window);


% Setup the text type for the window
% Screen('TextFont', window, 'Ariel');
Screen('TextSize', window, 36);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Here we set the size of the arms of our fixation cross
fixCrossDimPix = 40;


% Now we set the coordinates (these are all relative to zero we will let
% the drawing routine center the cross in the center of our monitor for us)
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];

% Set the line width for our fixation cross
lineWidthPix = 4;


%----------------------------------------------------------------------
%                       Keyboard information
%----------------------------------------------------------------------

% Define the keyboard keys that are listened for. We will be using the left
% and right arrow keys as response keys for the task and the escape key as
% a exit/reset key
escapeKey = KbName('ESCAPE');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');
downKey = KbName('DownArrow');
HideCursor;

% for the UG datafile
% num_errors is the total number of times they get one of the test
% questions incorrect, and total_time is the total amount of time they
% spent during the instructions phase.
sub_data          = NaN(1, 2);
sub_data_colnames = {'num_errors', 'total_time'};


%----------------------------------------------------------------------
%                       The instructions begin
%----------------------------------------------------------------------
txt        = sprintf('Please read the following very carefully, and ask the researcher if you have any questions.');
finished   = 0;
page_count = 1;
while finished == 0
    
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck;
    if keyIsDown == 1
        pause(.2)
    end
    
    if keyCode(leftKey) == 1 && page_count > 1
        page_count = page_count-1;
    elseif keyCode(rightKey) == 1
        page_count = page_count+1;
    end
    
    if page_count <= 1
        txt = sprintf(['Please read the following very carefully,\n'...
            'and ask the researcher if you have any questions.']);
    elseif page_count == 2
        txt = sprintf(['In this experiment, you will be completing 4 different tasks:\n'...
            'The ultimatum task, a task to assess your confidence in your\n'...
            'performance in the ultimatum task, an exchange task, and a\n'...
            'questionnaire. Instructions for each will be provided.\n\n']);
    elseif page_count == 3
        txt = sprintf(['One random trial from each of these tasks\n'...
            'will be chosen to calculate the extra money you earn,\n'...
            'therefore you will be paid out the extra amount earned on 4 trials.\n'...
            'This means that the decisions you make on each task will affect\n'...
            'your extra payout.']);
    elseif page_count == 4
        txt = sprintf('ULTIMATUM TASK');
    elseif page_count == 5
        txt = sprintf(['In the ultimatum task , you will be acting as\n'...
            'the so-called proposer. This means that, on every trial,\n'...
            'you will have a 20 monetary unit (MU) endowment (20 MU = %s2.50).\n'...
            'With this endowment, you will make a take-it-or-leave-it offer\n'...
            'to the other player - the responder.'], char(8364));
    elseif page_count == 6
        txt = sprintf(['The responder will either accept your offer,\n'...
            'in which case both you and the responder keep whatever\n'...
            'money you have proposed, or reject the offer, in which case\n'...
            'both you and the responder receive 0 MU for that trial.\n'...
            'Both your identity, and that of the responder,\n'...
            'will remain hidden']);
    elseif page_count == 7
        txt = sprintf(['As the proposer you can offer any amount\n'...
            'between 0 and 20 as proposal on how to divide the 20 MU''s\n'...
            'between you and the responder.']);
    elseif page_count == 8
        txt = sprintf([]);
        
    elseif page_count == 10
        txt = sprintf(['End of instructions.\n'...
            'Once you have moved past this page,\n'...
            'you will NOT be able to return. If anything is unclear,\n'...
            'please re-read the instructions now. Otherwise, continue.']);
    elseif page_count > 10
        finished = 1;
    end
    
    DrawFormattedText(window, txt, 'center', 'center');
    DrawFormattedText(window, 'press right key to continute, left key to go back', 'center', yCenter+windowRect(4)*.4);
    Screen('Flip', window, grey);
    %     KbStrokeWait;
    
end
sca;

%
%
%     txt = sprintf(['In this experiment, you will be completing 4 different tasks:\n'...
%         'The ultimatum task, a task to assess your confidence in your\n'...
%         'performance in the ultimatum task, an exchange task, and a\n'...
%         'questionnaire. Instructions for each will be provided.\n\n']);
%
%     DrawFormattedText(window, txt, 'center', 'center');
%     DrawFormattedText(window, 'press any key to continute', 'center', yCenter+windowRect(4)*.4);
%     Screen('Flip', window, grey);
%     KbStrokeWait;
%
%
%     txt = sprintf(['One random trial from each of these tasks\n'...
%         'will be chosen to calculate the extra money you earn,\n'...
%         'therefore you will be paid out the extra amount earned on 4 trials.\n'...
%         'This means that the decisions you make on each task will affect\n'...
%         'your extra payout.']);
%
%     DrawFormattedText(window, txt, 'center', 'center');
%     DrawFormattedText(window, 'press any key to continute', 'center', yCenter+windowRect(4)*.4);
%     Screen('Flip', window, grey);
%     KbStrokeWait;
