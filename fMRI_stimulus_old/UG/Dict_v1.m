% Dictator game
% Author: Michael Giffin, January, 2018
close all
% clearvars;

sub_num = input('What is the subject number?\n');

sca;

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
escapeKey   = KbName('ESCAPE');
leftKey     = KbName('LeftArrow');
rightKey    = KbName('RightArrow');
downKey     = KbName('DownArrow');

option1_key = KbName('a');
option2_key = KbName('l');
HideCursor;


%----------------------------------------------------------------------
%                       Number of trials and blocks
%----------------------------------------------------------------------
% first create all possible combinations of offers and shuffle them up real
% good
offer_choice = repmat(nchoosek(0:20,2),1,1);
% offer_choice = Shuffle(offer_choice);
offer_choice(1:length(offer_choice)/2, [1 2]) = offer_choice(1:length(offer_choice)/2, [2 1]);
offer_choice = Shuffle(offer_choice, 2);

% numTrials  = length(offer_choice); 
numTrials  = 9; 

%----------------------------------------------------------------------
%                       Pre-allocating data storage variables
%----------------------------------------------------------------------

% for the UG datafile
% column 1: offer
% column 2: option1
% column 3: option2
% column 4: RT
% % sub_data          = NaN(numTrials*numBlocks, 7);
sub_data          = NaN(numTrials, 4);
sub_data_colnames = {'offer', 'option1', 'option2', 'RT'};


instruction_txt1 = sprintf(['In this next section\n'...
    'you will be playing the exchange task.\n',...
    'In this task, similar to the ultimatum task,\n',...
    'you will be making offers to another individual.\n'...
    'However, in this task the other individual will have to\n',...
    'accept your offer. \n\nTherefore, whatever offer you make in\n',...
    'this task will be accepted.']);


instruction_txt2 = sprintf(['This task does not have a non-social condition,\n',...
    'each decision you make will affect another individual.\n'...
    'Importantly, during this task the other\n',...
    'individual whose payment depends on your offers\n'...
    'is not one of the individuals from the previous tasks.\n',...
    'The individuals you decide for in this task are\n'...
    'new individuals, all of whom will receive only what\n',...
    'you offer them (i.e. they have no extra money)']);



instruction_txt3 = sprintf(['You will be presented with a series of\n'...
    'binary choices from which you can select your offer.\n'...
    'For example, 10 and 2. If you select 10, then both\n'...
    'and the other player will receive 10. However,\n'...
    'if you select 2, the other player receives 2, and\n'...
    'and you receive 18']);


vbl = Screen('Flip', window);

% Length of time and number of frames we will use for each drawing test
numSecs = 1;
numFrames = round(numSecs / ifi);


DrawFormattedText(window, instruction_txt1, 'center', 'center');
DrawFormattedText(window, 'press any key to continute', 'center', yCenter+windowRect(4)*.4);
Screen('Flip', window, grey);
KbStrokeWait;

DrawFormattedText(window, instruction_txt2, 'center', 'center');
DrawFormattedText(window, 'press any key to continute', 'center', yCenter+windowRect(4)*.4);
Screen('Flip', window, grey);
KbStrokeWait;

DrawFormattedText(window, instruction_txt3, 'center', 'center');
DrawFormattedText(window, 'press any key to continute', 'center', yCenter+windowRect(4)*.4);
Screen('Flip', window, grey);
KbStrokeWait;

for trial = 1:numTrials
    
    % for getting RT
    T0 = GetSecs;
    
    answer = 0;
    while answer == 0
        
        % put the text over the image
        DrawFormattedText(window, 'Your offer', 'center', yCenter-windowRect(4)*.05);
        
        DrawFormattedText(window, [num2str(offer_choice(trial, 1))], xCenter-windowRect(3)*.2, yCenter+windowRect(4)*.05);
        DrawFormattedText(window, [num2str(offer_choice(trial, 2))], xCenter+windowRect(3)*.2, yCenter+windowRect(4)*.05);
        
        % Flip to the screen
        Screen('Flip', window, grey);
        
        [keyIsDown, secs, keyCode, deltaSecs] = KbCheck;
        
        if keyCode(option1_key) == 1
            answer = 1;
            offer = offer_choice(trial, 1);
            
            RT = secs - T0;
            for frame = 1:numFrames/2
                Screen('TextStyle', window, 1);
                
                DrawFormattedText(window, [num2str(offer_choice(trial, 1))], xCenter-windowRect(3)*.2, yCenter+windowRect(4)*.05);
                
                Screen('Flip', window, grey);
                Screen('TextStyle', window, 0);
            end
        elseif keyCode(option2_key) == 1
            answer = 1;
            offer = offer_choice(trial, 2);
            
            RT = secs - T0;
            for frame = 1:numFrames/2
                Screen('TextStyle', window, 1);
                
                DrawFormattedText(window, [num2str(offer_choice(trial, 2))], xCenter+windowRect(3)*.2, yCenter+windowRect(4)*.05);
                
                Screen('Flip', window, grey);
                Screen('TextStyle', window, 0);
            end
        end
    end
    sub_data(trial, 1) = offer;
    sub_data(trial, 2) = offer_choice(trial, 1);
    sub_data(trial, 3) = offer_choice(trial, 2);
    sub_data(trial, 4) = RT;
end

% first check to see if there is already a file for this subject and this
% block and if so, save the file with the date at the end
if sub_num > 99 
    if exist(sprintf('data%spilot%ssub%d_dict.mat', filesep, filesep, sub_num), 'file') == 2
        save(sprintf('data%spilot%ssub%d_dict_%s.mat', filesep, filesep, sub_num, date), 'sub_data', 'sub_data_colnames');
    else
        save(sprintf('data%spilot%ssub%d_dict.mat', filesep, filesep, sub_num), 'sub_data', 'sub_data_colnames');
    end
elseif sub_num > 9 && sub_num < 100
    if exist(sprintf('data%spilot%ssub0%d_dict.mat', filesep, filesep, sub_num), 'file') == 2
        save(sprintf('data%spilot%ssub0%d_dict_%s.mat', filesep, filesep, sub_num, date), 'sub_data', 'sub_data_colnames');
    else
        save(sprintf('data%spilot%ssub0%d_dict.mat', filesep, filesep, sub_num), 'sub_data', 'sub_data_colnames');
    end
elseif sub_num < 10
    if exist(sprintf('data%spilot%ssub00%d_dict.mat', filesep, filesep, sub_num), 'file') == 2
        save(sprintf('data%spilot%ssub00%d_dict_%s.mat', filesep, filesep, sub_num, date), 'sub_data', 'sub_data_colnames');
    else
        save(sprintf('data%spilot%ssub00%d_dict.mat', filesep, filesep, sub_num), 'sub_data', 'sub_data_colnames');
    end
end

sca;