% Ultimatum Game for social learning in fMRI
% this particular version is meant to be run during the acquisition of the
% T1, meaning we won't need to collect any concurrent neural data with this
% tastk. And it's 
% Author: Michael Giffin, January, 2018
%------ inputs -------

% stuff the script will ask me for before running the experiment
sub_num = input('What is the subject number?\n');
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
condition = input('Condition (1:4): \n');
% condition = 1;

% this will be given to the subject twice, once after the first two blocks,
% and once after the second two blocks
block = input('Block (1:2): \n');
% block = 1;
%------ inputs -------

% Clear the workspace and the screen
sca;
close all;
% clearvars;

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


%=========================================================================
%               OF CRITICAL IMPORTANCE
%               wait here for the first pulse 
%               of the fMRI signal
%=========================================================================
% Wait for the fMRI pulse, need to know what this is registered as, can
% either be a serial port or registered as a particulary key. 
% In Leiden it is registered as the number 5
% IMPORTANTLY, this is 5% on my computer
fMRI_sim = KbName('5%');
fMRI_key = KbName('5');
% fMRI_key = KbName('t');
fMRI_wait = 0;
while fMRI_wait == 0
    DrawFormattedText(window, 'Please wait...', 'center', 'center');
    Screen('Flip', window, grey);
    [keyIsDown, secs, keyCode, deltaSecs] = KbCheck;
    if keyCode(fMRI_key) == 1 || keyCode(fMRI_sim) == 1
        %     if keyCode(fMRI_key) == 1
        fMRI_wait = 1;
        T0 = 0;
        tic
    end
end


%----------------------------------------------------------------------
%                       Timing Information
%----------------------------------------------------------------------

% Interstimulus interval time in seconds and frames
isiTimeSecs = 1;
isiTimeFrames = round(isiTimeSecs / ifi);

% Numer of frames to wait before re-drawing
% Here we use to a waitframes number greater then 1 to flip at a rate not
% equal to the monitors refreash rate. For this example, once per second,
% to the nearest frame
flipSecs = 10;
waitframes = round(flipSecs / ifi);



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


%----------------------------------------------------------------------
%                       Number of trials and blocks
%----------------------------------------------------------------------

% numTrials  = 72; % must be divisible by 3, since that's how many opponents are in each block
numTrials  = 3; % must be divisible by 3, since that's how many opponents are in each block
%----------------------------------------------------------------------
%                       Pre-allocating data storage variables
%----------------------------------------------------------------------

% for the UG datafile
% column 1: offer
% column 2: accepted (1 for yes, 0 for no)
% column 3: payment
% column 4: social (1 for social, 0 for non-social)
% column 5: opponent, 0 for starting endowment of 0, 1 for starting
%           endowment of 10, and 2 for starting endowment of 20
% column 6: RT
% column 7: fMRI pulse?
% % sub_data          = NaN(numTrials*numBlocks, 7);
sub_data          = NaN(numTrials, 7);
sub_data_colnames = {'offer', 'accepted', 'payment', 'social', 'opponent', 'RT', 'offer_start_position'};

% cell array to store shapes that represented each opponent
% column 1: shape
% column 2: opponent that the shape represents
% % sub_opp_shapes  = cell(numTrials*numBlocks, 2);
sub_opp_shapes  = cell(numTrials, 2);

% tells me the onset of each window in seconds since the first pulse
% recieved from the scanner
window_times = struct('fixation', NaN(numTrials, 1), 'wait1', NaN(numTrials, 1),...
    'opponent', NaN(numTrials, 1), 'wait2', NaN(numTrials, 1),...
    'offer_start', NaN(numTrials, 1), 'offer_end', NaN(numTrials, 1),...
    'wait3_post_offer', NaN(numTrials, 1), 'accept', NaN(numTrials, 1),...
    'wait4', NaN(numTrials, 1), 'payment', NaN(numTrials, 1),...
    'wait5', NaN(numTrials, 1));

%----------------------------------------------------------------------
%                   stuff for the slider
%----------------------------------------------------------------------
% bgColor   = [128 128 128];
% [myScreen, rect] = Screen('OpenWindow', screenNumber, bgColor);
% 
% center           = round([rect(3) rect(4)]/2);
% 
question  = 'Probability of acceptance';
endPoints = {'no', 'yes'};


%----------------------------------------------------------------------
%           Fixation cross stuff
%----------------------------------------------------------------------

% Draw the fixation cross in white, set it to the center of our screen and
% set good quality antialiasing
Screen('DrawLines', window, allCoords,...
    lineWidthPix, white, [xCenter yCenter], 2);


vbl = Screen('Flip', window);

% Length of time and number of frames we will use for each drawing test
numSecs = 1;
numFrames = round(numSecs / ifi);

% some stuff for the scalar and it's relation to the size of the screen in
% terms of ratio
scalaPosition = .8;
scalaLength   = .9;
%==========================================================================
%----------------------------------------------------------------------
%          Opponents acceptance functions and shuffled order
%----------------------------------------------------------------------

% setting up opponents 
endow0   =  [0.0707, 0.1032, 0.1482, 0.2083, 0.2846, 0.3757, 0.4764, 0.5791, 0.6754, 0.7589, 0.8264, 0.8780, 0.9159, 0.9427, 0.9614, 0.9741, 0.9827, 0.9885, 0.9924, 0.9949, 0.9967];
endow10  =  [0.2844, 0.3839, 0.4942, 0.6051, 0.7061, 0.7902, 0.8552, 0.9026, 0.9356, 0.9579, 0.9728, 0.9825, 0.9887, 0.9928, 0.9954, 0.9971, 0.9981, 0.9988, 0.9992, 0.9995, 0.9997];
endow20  =  [0.8268, 0.8855, 0.9261, 0.9530, 0.9705, 0.9816, 0.9885, 0.9929, 0.9956, 0.9973, 0.9983, 0.9990, 0.9994, 0.9996, 0.9998, 0.9998, 0.9999, 0.9999, 1.0000, 1.0000, 1.0000];
% % opp_txt  = 'You are playing against an opponent from';

opponents   = [zeros(numTrials/3, 1); ones(numTrials/3, 1); repmat(2, [numTrials/3, 1])]; % numTrials/3 because I have 3 opponents per session
opponents   = Shuffle(opponents);

%----------------------------------------------------------------------
%         assigning conditions 
%----------------------------------------------------------------------
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
elseif condition == 4
    socmat = [0 1 0 1];
    socpo  = 'frame';
    nonsocpo= 'filled';
end


%----------------------------------------------------------------------
%          Opponents images
%----------------------------------------------------------------------
% opponent images
if block == 1
    block = 1:2;
elseif block == 2
    block = 3:4;
else
    sca;
    error('Blocks need to be either 1 or 2 for this task');
end


% % imstruct        = dir(['block', num2str(block), '_images' filesep '*png*']);

% the image directories are set up such that each shape appears twice, once
% as a frame and once filled. Each shape will be assigned to represent one
% opponent, and will never be repeated. So if a framed triangle is selected
% to represent a nonsocial opponent, a filled triangle will never be used
% for a social opponent


all_count = 0; % the counter that keeps trac so I can record to behavioral data
               % was originally here because I had a loop over blocks as
               % well as trials, but now each block is saved independently 
%==========================================================================
%                   The experiment loop
%==========================================================================
opponents   = Shuffle(opponents); % shuffle the opponents on each block so subjects never play against the same order
% % 
% % social = socmat(block);
% % if social == 1
% %     % %         imstruct = soc_img_dir;
% %     soctxt   = 'social';
% % else
% %     % %         imstruct = nonsoc_img_dir;
% %     soctxt   = 'non-social';
% % end

instruction_txt = sprintf(['In the following section\n'...
    'you will be asked about how likely you think it is\n'...
    'that a particular opponent against whom you just played\n'...
    'will accept or reject your offer\n\n'...
    'Please let the experimenter know\n'...
    'If you have any questions about this task\n\n'...
    'Otherwise, press any key to continue']);

DrawFormattedText(window, instruction_txt, 'center', 'center');
Screen('Flip', window, grey);

KbStrokeWait;

for all_blocks = 1:length(block)
    
    imstruct        = dir(['block', num2str(block(all_blocks)), '_images' filesep '*png*']);
    
    for trial = 1:numTrials
        all_count = all_count + 1;
        % selecting the opponent against which the subject is playing
        if opponents(trial) == 0
            [img, ~, alpha] = imread(['block', num2str(block(all_blocks)), '_images', filesep, imstruct(1).name]);
            img_name = imstruct(1).name;
        elseif opponents(trial) == 1
            [img, ~, alpha] = imread(['block', num2str(block(all_blocks)), '_images', filesep, imstruct(2).name]);
            img_name = imstruct(2).name;
        elseif opponents(trial) == 2
            [img, ~, alpha] = imread(['block', num2str(block(all_blocks)), '_images', filesep, imstruct(3).name]);
            img_name = imstruct(3).name;
        end
        
        img = double(img);
        img = rgb2gray(img);
        img(img == 1) = grey;
        % img(:, :, 4) = alpha;
        texture2 = Screen('MakeTexture', window, img);
        
        for all_offers = 1:21 % ask about each offer for each opponent
            
% %             window_times.opponent(trial) = toc;
% %             % screen indicating which opponent subject is playing against, ~2 seconds
% %             for frame = 1:numFrames+(numFrames/2)*rand
% %                 
% %                 % this puts it right in the middle of the screen
% %                 Screen('DrawTexture', window, texture2, [], [xCenter-windowRect(3)*.2, yCenter-windowRect(4)*.2, xCenter+windowRect(3)*.2, yCenter+windowRect(4)*.2]);
% %                 
% %                 % put the text over the image
% %                 DrawFormattedText(window, 'opponent from', 'center', yCenter-windowRect(4)*.17);
% %                 
% %                 DrawFormattedText(window, [soctxt '\noffers ' num2str(all_offers)], 'center', yCenter+windowRect(4)*.2);
% %                 
% %                 % Flip to the screen
% %                 Screen('Flip', window, grey);
% %                 
% %             end
            
            window_times.offer_start(trial) = toc;
            % the slider, where the decision is made
%             [offer, RT, answer, offer_start_position] = slideScale_probabilites(all_offers, social, window, ['offer of ', num2str(all_offers-1), '\n\nProbability of acceptance (in %)'], windowRect, endPoints, 'device', 'keyboard', 'scalaposition', scalaPosition, 'startposition', 'right', 'displayposition', false, 'image', img);
            [offer, RT, answer, offer_start_position] = slideScale_probabilites(all_offers, social, window, ['\n\nProbability of acceptance (in %)'], windowRect, endPoints, 'device', 'keyboard', 'scalaposition', scalaPosition, 'startposition', 'right', 'displayposition', false, 'image', img);
            
        end
    end
end



% % % from full UG
% % for trial = 1:numTrials
% %     all_count = all_count + 1;
% %     
% %     if trial == 1
% %         DrawFormattedText(window, instruction_txt, 'center', 'center');
% %         Screen('Flip', window, grey);
% %         
% %         KbStrokeWait;
% %     end
% %     
% % % %     window_times.fixation(trial) = toc;
% % % %     % fixation cross loop, jittered between .5 and 1.5 seconds
% % % %     for frame = 1:(numFrames/2) + numFrames*rand
% % % %         % Draw the fixation cross in white, set it to the center of our screen and
% % % %         % set good quality antialiasing
% % % %         Screen('DrawLines', window, allCoords,lineWidthPix, white, [xCenter yCenter], 2);
% % % %         
% % % %         % Flip to the screen
% % % %         Screen('Flip', window);
% % % %         
% % % %     end
% % % %     
% % % %     window_times.wait1(trial) = toc;
% % % %     % blank wait screen between .5 and 1 seconds
% % % %     for frame = 1:(numFrames/2)+(numFrames/2)*rand
% % % %         Screen('Flip', window, grey);
% % % %     end
% %     
% %     % selecting the opponent against which the subject is playing
% %     if opponents(trial) == 0
% %         [img, ~, alpha] = imread(['block', num2str(block), '_images', filesep, imstruct(1).name]);
% %         img_name = imstruct(1).name;
% %     elseif opponents(trial) == 1
% %         [img, ~, alpha] = imread(['block', num2str(block), '_images', filesep, imstruct(2).name]);
% %         img_name = imstruct(2).name;
% %     elseif opponents(trial) == 2
% %         [img, ~, alpha] = imread(['block', num2str(block), '_images', filesep, imstruct(3).name]);
% %         img_name = imstruct(3).name;
% %     end
% %     
% %     img = double(img);
% %     img = rgb2gray(img);
% %     img(img == 1) = grey;
% %     % img(:, :, 4) = alpha;
% %     texture2 = Screen('MakeTexture', window, img);
% %     
% %     
% %     window_times.opponent(trial) = toc;
% %     % screen indicating which opponent subject is playing against, ~2 seconds
% %     for frame = 1:numFrames+(numFrames/2)*rand
% %             
% %         % this puts it right in the middle of the screen
% %         Screen('DrawTexture', window, texture2, [], [xCenter-windowRect(3)*.2, yCenter-windowRect(4)*.2, xCenter+windowRect(3)*.2, yCenter+windowRect(4)*.2]);
% %         
% %         % put the text over the image
% %         DrawFormattedText(window, 'opponent from', 'center', yCenter-windowRect(4)*.17);
% %         
% %         DrawFormattedText(window, soctxt, 'center', yCenter+windowRect(4)*.2);
% %         
% %         % Flip to the screen
% %         Screen('Flip', window, grey);
% %         
% %     end
% %     
% %     window_times.wait2(trial) = toc;
% %     % blank wait screen between .5 and 1 seconds
% %     for frame = 1:(numFrames/2)+(numFrames/2)*rand
% %         Screen('Flip', window, grey);
% %     end
% %     
% %     window_times.offer_start(trial) = toc;
% %     % the slider, where the decision is made
% %     [offer, RT, answer, offer_start_position] = slideScale_probabilites(window, question, windowRect, endPoints, 'device', 'keyboard', 'scalaposition', scalaPosition, 'startposition', 'right', 'displayposition', false);
% %     window_times.offer_end(trial) = toc;
% %     time_since_pulse              = window_times.offer_end(trial);
% %     
% %     window_times.wait3_post_offer(trial) = toc;
% %     % blank wait screen between 1.5 and 2.5 seconds
% %     for frame = 1:(numFrames/2)+(numFrames/2)+(numFrames/2)*rand
% %         Screen('Flip', window, grey);
% %         
% %     end
% %     
% %     % calculate whether or not the offer was accepted by the opponent
% %     if opponents(trial) == 0
% %         accept = rand < endow0(offer+1); % +1 because I can't index with 0
% %     elseif opponents(trial) == 1
% %         accept = rand < endow10(offer+1); % +1 because I can't index with 0
% %     elseif opponents(trial) == 2
% %         accept = rand < endow20(offer+1); % +1 because I can't index with 0
% %     end
% %     
% %     
% %     if accept == 1
% %         acc_txt   = sprintf('ACCEPTED');
% %         acc_color = [0 255 0]; % green
% %     else
% %         acc_txt = sprintf('REJECTED');
% %         acc_color = [255 0 0]; % red
% %     end
% %     
% %     
% %     window_times.accept(trial) = toc;
% %     % screen indicating whether the offer was accepted or not, 2
% %     % seconds + jitter around .5 seconds
% %     for frame = 1:numFrames+(numFrames/2)*rand
% %         % this puts it right in the middle of the screen
% %         Screen('DrawTexture', window, texture2, [], [xCenter-windowRect(3)*.2, yCenter-windowRect(4)*.2, xCenter+windowRect(3)*.2, yCenter+windowRect(4)*.2]);
% %         
% %         DrawFormattedText(window, acc_txt, 'center', yCenter+windowRect(4)*.2, acc_color);
% %         DrawFormattedText(window, 'Opponent', 'center',yCenter-windowRect(4)*.17, [0 0 0]);
% %         
% %         % Flip to the screen
% %         Screen('Flip', window, grey);
% %     end
% %     
% % % %     window_times.wait4(trial) = toc;
% % % %     % blank wait screen between .5 and 1 seconds
% % % %     for frame = 1:(numFrames/4)+(numFrames/2)*rand
% % % %         Screen('Flip', window, grey);
% % % %     end
% %     
% %     if accept == 1
% %         payment = 20-offer;
% %         acc_txt = sprintf('%d', payment);
% % %         acc_txt = sprintf('\n\n%d', payment);
% % % %         acc_color = [0 255 0]; % green
% %     else
% %         payment = 0;
% %         acc_txt = sprintf('%d', payment);
% % %         acc_txt = sprintf('\n\n%d', payment);
% % % %         acc_color = [255 0 0]; % red
% %     end
% %     
% %     window_times.payment(trial) = toc;
% %     % screen indicating how much the subject received
% %     % still need to tweak it to make the amount they earned bigger
% %     for frame = 1:(numFrames+(numFrames/2))+(numFrames/4)*rand
% %         
% % %         DrawFormattedText(window, ['\n\n', acc_txt], 'center','center');
% % %         DrawFormattedText(window, 'You gain', 'center','center');
% % %         
% %         
% %         % this puts it right in the middle of the screen
% %         Screen('DrawTexture', window, texture2, [], [xCenter-windowRect(3)*.2, yCenter-windowRect(4)*.2, xCenter+windowRect(3)*.2, yCenter+windowRect(4)*.2]);
% %         
% %         DrawFormattedText(window, ['You receive ', acc_txt], 'center', yCenter+windowRect(4)*.2);
% %         DrawFormattedText(window, ['You offered ', num2str(offer)], 'center',yCenter-windowRect(4)*.17);
% %         
% %         Screen('Flip', window, grey);
% %     end
% %     
% %     window_times.wait5(trial) = toc;
% %     % blank wait screen between .5 and 1 seconds
% %     for frame = 1:(numFrames/2)+(numFrames/2)*rand
% %         Screen('Flip', window, grey);
% %     end
% %     
% %     sub_data(all_count, 1)       = offer;
% %     sub_data(all_count, 2)       = accept;
% %     sub_data(all_count, 3)       = payment;
% %     sub_data(all_count, 4)       = social;
% %     sub_data(all_count, 5)       = opponents(trial);
% %     sub_data(all_count, 6)       = RT;
% %     sub_data(all_count, 7)       = offer_start_position;
% %     
% %     sub_opp_shapes{all_count, 1} = img_name;
% %     sub_opp_shapes{all_count, 2} = opponents(trial);
% %     
% % end

% for the data files that I'm saving
if social == 1
    soctxt   = 'social';
else
    soctxt   = 'nonsoc';
end

% first check to see if there is already a file for this subject and this
% block and if so, save the file with the date at the end
if sub_num > 99 
    if exist(sprintf('data%spilot%ssub_prob%d_%s_%d.mat', filesep, filesep, sub_num, soctxt, block), 'file') == 2
        save(sprintf('data%spilot%ssub_prob%d_%s_%d_%s.mat', filesep, filesep, sub_num, soctxt, block, date), 'sub_data', 'condition', 'sub_opp_shapes', 'sub_data_colnames', 'window_times', 'soctxt', 'block');
    else
        save(sprintf('data%spilot%ssub_prob%d_%s_%d.mat', filesep, filesep, sub_num, soctxt, block), 'sub_data', 'condition', 'sub_opp_shapes', 'sub_data_colnames', 'window_times', 'soctxt', 'block');
    end
elseif sub_num > 9 && sub_num < 100
    if exist(sprintf('data%spilot%ssub0_prob%d_%s_%d.mat', filesep, filesep, sub_num, soctxt, block), 'file') == 2
        save(sprintf('data%spilot%ssub0_prob%d_%s_%d_%s.mat', filesep, filesep, sub_num, soctxt, block, date), 'sub_data', 'condition', 'sub_opp_shapes', 'sub_data_colnames', 'window_times', 'soctxt', 'block');
    else
        save(sprintf('data%spilot%ssub0_prob%d_%s_%d.mat', filesep, filesep, sub_num, soctxt, block), 'sub_data', 'condition', 'sub_opp_shapes', 'sub_data_colnames', 'window_times', 'soctxt', 'block');
    end
elseif sub_num < 10
    if exist(sprintf('data%spilot%ssub00_prob%d_%s_%d.mat', filesep, filesep, sub_num, soctxt, block), 'file') == 2
        save(sprintf('data%spilot%ssub00_prob%d_%s_%d_%s.mat', filesep, filesep, sub_num, soctxt, block, date), 'sub_data', 'condition', 'sub_opp_shapes', 'sub_data_colnames', 'window_times', 'soctxt', 'block');
    else
        save(sprintf('data%spilot%ssub00_prob%d_%s_%d.mat', filesep, filesep, sub_num, soctxt, block), 'sub_data', 'condition', 'sub_opp_shapes', 'sub_data_colnames', 'window_times', 'soctxt', 'block');
    end
end

% Wait for a key press
% KbStrokeWait;

% Clear the screen
sca;