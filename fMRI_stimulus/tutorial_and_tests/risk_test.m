
% just for testing out my slideScale_risk function
close all


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


% Length of time and number of frames we will use for each drawing test
numSecs = 1;
numFrames = round(numSecs / ifi);


%---------------%---------------%---------------%---------------%---------------
% the risk questionnaire
instruction_1 = sprintf(['Now you will be answer the risk questionnaire.\n',...
    'You have to choose one and only one of the following six gambles.\n'...
    'Each gamble has two possible outcomes (Roll High or Roll Low),\n'...
    'with the indicated probabilities of occurring. Your payment will\n'...
    'be determined by which of the options you select, and which of\m'...
    'the outcomes occurs.\n']);


DrawFormattedText(window, instruction_1, 'center', 'center');
DrawFormattedText(window, 'press any key to continute', 'center', yCenter+windowRect(4)*.4);
Screen('Flip', window, grey);
KbStrokeWait;

[img, ~, alpha] = imread(['~/Dropbox/RL/fMRI_stimulus/UG/risk_table' filesep 'risk_table.png']);

[risk_gamble, RT, answer, offer_start_position] = slideScale_risk(numFrames, window, question, windowRect, endPoints, 'device', 'keyboard', 'scalaposition', scalaPosition, 'startposition', 'right', 'displayposition', false, 'image', img);
