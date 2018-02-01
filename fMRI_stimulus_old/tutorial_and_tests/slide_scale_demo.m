% Preliminary stuff
% Clear Matlab/Octave window:
clc;


% check for Opengl compatibility, abort otherwise:
AssertOpenGL;

bgColor   = [128 128 128];

% HideCursor;

% Get information about the screen and set general things
Screen('Preference', 'SuppressAllWarnings',0);
Screen('Preference', 'SkipSyncTests', 1);

rect          = Screen('Rect',0);
screenRatio   = rect(3)/rect(4);
pixelSizes    = Screen('PixelSizes', 0);
startPosition = round([rect(3)/2, rect(4)/2]);


%----------------------------------------------
% Attempt to get screen to external secondary monitor
screenNumber = max(Screen('Screens'));

% Define black, white and grey
white = WhiteIndex(screenNumber);
grey = white / 2;
black = BlackIndex(screenNumber);

% Open the screen
% [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [], 32, 2);

%----------------------------------------------

% Creating screen etc.
% [myScreen, rect] = Screen('OpenWindow', 0, bgColor);
% draws to external screen
[myScreen, rect] = Screen('OpenWindow', screenNumber, bgColor);


center           = round([rect(3) rect(4)]/2);


question  = 'Select investment';
endPoints = {'no', 'yes'};

[position, RT, answer] = slideScale(myScreen, question, rect, endPoints, 'device', 'keyboard', 'scalaposition', 0.5, 'startposition', 'right', 'displayposition', true);


Screen('CloseAll')