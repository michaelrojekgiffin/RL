function [position, RT, answer, offer_start_position] = slideScale(numFrames, screenPointer, question, rect, endPoints, varargin)
%SLIDESCALE This funtion draws a slide scale on a PSYCHTOOLOX 3 screen and returns the
% position of the slider spaced between -100 and 100 as well as the rection time and if an answer was given.
%
%   Usage: [position, secs] = slideScale(ScreenPointer, question, center, rect, endPoints, varargin)
%   Mandatory input:
%    ScreenPointer  -> Pointer to the window.
%    question       -> Text string containing the question.
%    rect           -> Double contatining the screen size.
%                      Obtained with [myScreen, rect] = Screen('OpenWindow', 0);
%    endPoints      -> Cell containg the two text string of the left and right
%                      end of the scala. Exampe: endPoints = {'left, 'right'};
%
%   Varargin:
%    'linelength'     -> An integer specifying the lengths of the ticks in
%                        pixels. The default is 10.
%    'width'          -> An integer specifying the width of the scala line in
%                        pixels. The default is 3.
%    'startposition'  -> Choose 'right', 'left' or 'center' start position.
%                        Default is center.
%    'scalalength'    -> Double value between 0 and 1 for the length of the
%                        scale. The default is 0.9.
%    'scalaposition'  -> Double value between 0 and 1 for the position of the
%                        scale. 0 is top and 1 is bottom. Default is 0.8.
%    'device'         -> A string specifying the response device. Either 'mouse' 
%                        or 'keyboard'. The default is 'mouse'.
%    'responsekey'    -> String containing name of the key from the keyboard to log the
%                        response. The default is 'return'.
%    'slidecolor'     -> Vector for the color value of the slider [r g b] 
%                        from 0 to 255. The dedult is red [255 0 0].
%    'scalacolor'     -> Vector for the color value of the scale [r g b] 
%                        from 0 to 255.The dedult is black [0 0 0].
%    'aborttime'      -> Double specifying the time in seconds after which
%                        the function should be aborted. In this case no
%                        answer is saved. The default is 8 secs.
%    'image'          -> An image saved in a uint8 matrix. Use
%                        imread('image.png') to load an image file.
%    'displaypoition' -> If true, the position of the slider is displayed. 
%                        The default is false. 
%
%   Output:
%    'position'      -> Deviation from zero in percentage, 
%                       with -100 <= position <= 100 to indicate left-sided
%                       and right-sided deviation.
%    'RT'            -> Reaction time in milliseconds.
%    'answer'        -> If 0, no answer has been given. Otherwise this
%                       variable is 1.
%
%   Original Author: Joern Alexander Quent
%   Tinkerer       : Michael Giffin (made specific for keyboard movement and
%                                    multiple ticks on scalar for UG)
%   e-mail: alexander.quent@rub.de
%   e-mail: m.r.giffin@fsw.leidenuniv.nl
%   Version history:
%                    1.0 - 4. January 2016 - First draft
%                    1.1 - 18. Feburary 2016 - Added abort time and option to
%                    choose between mouse and key board
%                    1.2 - 5. October 2016 - End points will be aligned to end
%                    ticks
%                    1.3 - 06/01/2017 - Added the possibility to display an
%                    image
%                    1.4 - 5. May 2017 - Added the possibility to choose a
%                    start position
%                    1.5 - 7. November 2017 - Added the possibility to display
%                    the position of the slider under the scale.
%                    1.6 - 27. November 2017 - The function now waits until
%                    all keys are released before exiting. 
%                    1.7 - 28. November 2017 - More than one screen
%                    supported now.
%                    1.8 - 29. November 2017 - Fixed issue that mouse is
%                    not properly in windowed mode.
%                    1.9 - 7. December 2017 - If an image is drawn, the
%                    corresponding texture is deleted at the end.

% Michael Giffin enters the scene
%                    2.0 - 20. December 2017 - change scale movability from
%                    mouse to keyboard
%                    2.1 - 2. January 2018 - change scale to have 20 tick
%                    marks for offer selection for ultimatum game



%% Parse input arguments
% Default values
center        = round([rect(3) rect(4)]/2);
lineLength    = 10;
width         = 3;
scalaLength   = 0.9;
scalaPosition = 0.8;
sliderColor    = [255 0 0];
scaleColor    = [0 0 0];
device        = 'mouse';
aborttime     = 5;
% aborttime     = .5;
% responseKey   = KbName('space');

responseKey   = KbName('6^');

GetMouseIndices;
drawImage     = 0;
startPosition = 'center';
displayPos    = false;

i = 1;
while(i<=length(varargin))
    switch lower(varargin{i})
        case 'linelength'
            i             = i + 1;
            lineLength    = varargin{i};
            i             = i + 1;
        case 'width'
            i             = i + 1;
            width         = varargin{i};
            i             = i + 1;
        case 'startposition'
            i             = i + 1;
            startPosition = varargin{i};
            i             = i + 1;
        case 'scalalength'
            i             = i + 1;
            scalaLength   = varargin{i};
            i             = i + 1;
        case 'scalaposition'
            i             = i + 1;
            scalaPosition = varargin{i};
            i             = i + 1;
        case 'device' 
            i             = i + 1;
            device = varargin{i};
            i             = i + 1;
        case 'responsekey'
            i             = i + 1;
            responseKey   = KbName(varargin{i});
            i             = i + 1;
        case 'slidecolor'
            i             = i + 1;
            sliderColor    = varargin{i};
            i             = i + 1;
        case 'scalacolor'
            i             = i + 1;
            scaleColor    = varargin{i};
            i             = i + 1;
        case 'aborttime'
            i             = i + 1;
            aborttime     = varargin{i};
            i             = i + 1;
        case 'image'
            i             = i + 1;
            image         = varargin{i};
            i             = i + 1;
            imageSize     = size(image);
            stimuli       = Screen('MakeTexture', screenPointer, image);
            drawImage     = 1; 
        case 'displayposition'
            i             = i + 1;
            displayPos    = varargin{i};
            i             = i + 1;
    end
end

% Sets the default key depending on choosen device
if strcmp(device, 'mouse')
    responseKey   = 1; % X mouse button
end

%% Checking number of screens and parsing size of the global screen
screens       = Screen('Screens');
if length(screens) > 1 % Checks for the number of screens
    screenNum        = 1;
else
    screenNum        = 0;
end
globalRect          = Screen('Rect', screenNum);

%% Coordinates of scale lines and text bounds
if strcmp(startPosition, 'right')
    x = globalRect(3)*scalaLength;
    position = 20;
elseif strcmp(startPosition, 'center')
    x = globalRect(3)/2;
elseif strcmp(startPosition, 'left')
    x = globalRect(3)*(1-scalaLength);
    position = 0;
else
    error('Only right, center and left are possible start positions');
end


% SetMouse(round(x), round(rect(4)*scalaPosition), screenPointer, 1);
midTick    = [center(1), rect(4)*scalaPosition - lineLength - (lineLength/2), center(1), rect(4)*scalaPosition + lineLength + (lineLength/2)];
leftTick   = [rect(3)*(1-scalaLength), rect(4)*scalaPosition - lineLength, rect(3)*(1-scalaLength), rect(4)*scalaPosition  + lineLength];
rightTick  = [rect(3)*scalaLength, rect(4)*scalaPosition - lineLength rect(3)*scalaLength, rect(4)*scalaPosition  + lineLength];
horzLine   = [rect(3)*scalaLength, rect(4)*scalaPosition, rect(3)*(1-scalaLength), rect(4)*scalaPosition];
textBounds = [Screen('TextBounds', screenPointer, endPoints{1}); Screen('TextBounds', screenPointer, endPoints{2})];


if drawImage == 1
    rectImage  = [center(1) - imageSize(2)/2 rect(4)*(scalaPosition - 0.2) - imageSize(1) center(1) + imageSize(2)/2 rect(4)*(scalaPosition - 0.2)];
    if rect(4)*(scalaPosition - 0.2) - imageSize(1) < 0
        error('The height of the image is too large. Either lower your scale or use the smaller image.');
    end
end

% Calculate the range of the scale, which will be need to calculate the
% position
scaleRange        = round(rect(3)*(1-scalaLength)):round(rect(3)*scalaLength); % Calculates the range of the scale
scaleRangeShifted = round((scaleRange)-mean(scaleRange));                      % Shift the range of scale so it is symmetrical around zero


% my attempt to get ticks for all 21 offers
all_ticks = NaN(21, 4);
xticker = 0;
for tickem = 1:length(all_ticks)
    all_ticks(tickem, :)  = [rect(3)*(1-scalaLength) + xticker, rect(4)*scalaPosition - lineLength, rect(3)*(1-scalaLength) + xticker, rect(4)*scalaPosition  + lineLength];
    xticker = xticker + (length(scaleRange)/20);
end

% start in random point on scale for each trial within range in the middle
% quadrant of the scale
rnd_range = [ceil(length(all_ticks)/4), ceil(length(all_ticks) - (length(all_ticks)/4))];
x         = datasample(all_ticks(rnd_range(1):rnd_range(end), 1), 1);
position  = find(x(:, 1) == all_ticks(:, 1))-1;

offer_start_position = position;

SetMouse(round(x), round(rect(4)*scalaPosition), screenPointer, 1);


%% HILMAR
% F = figure;
% % text(-1,0,{'Glorious Experiment, she is now runnink, yes?', 'Please to be pressink "5", I will collectink data.', '', 'Please to be pressink key of escapeness when finishedski, da?'}, 'fontsize', 14);
% % axis off; axis ([-1, 1, -1, 1]);
% set (F, 'keypressfcn', @detectKeyPresses);
% waitfor(F);
% 
%  function detectKeyPresses (CallerHandle, KeyPressEvent, StartTime)
%     if strcmp (KeyPressEvent.Key, '1'); 
% %    if strcmp (KeyPressEvent.Key, '5'); 
%      global fMRIData; fMRIData(end+1) = toc(StartTime);
%      fprintf('"5" pressed at %d seconds.\n', fMRIData(end)); return
%     elseif strcmp (KeyPressEvent.Key, '2'); 
% %    if strcmp (KeyPressEvent.Key, '5'); 
%      global fMRIData; fMRIData(end+1) = toc(StartTime);
%      fprintf('"5" pressed at %d seconds.\n', fMRIData(end)); return
%    elseif strcmp (KeyPressEvent.Key, 'escape');
%      disp ('Escape Pressed. Ending Experiment'); 
%      close (CallerHandle);
%     end
%  end
%% Loop for scale loop
t0                         = GetSecs;
answer                     = 0;

% assigning the keyboard namesS
leftKey = KbName('1!');
rightKey = KbName('2@');
five = KbName('5');
fivefive = KbName('5%');
% 
% leftKey = KbName('LeftArrow');
% rightKey = KbName('RightArrow');
desperate = 0;

keyscroll_helper  = 0; % value to allow rapid scrolling
keyscroll_tracker = [];
while answer == 0
     % I think this is where I need to change things so that it's not the
    % mouse that's changing the position of the scaler, but the keys,
    % because eventually I'm going to have to change to the fMRI selectors.
    
%     %original
%     [x,y,buttons,focus,valuators,valinfo] = GetMouse(screenPointer, 1);
% %     the original
%     if x > rect(3)*scalaLength
%         x = rect(3)*scalaLength;
%     elseif x < rect(3)*(1-scalaLength)
%         x = rect(3)*(1-scalaLength);
%     end
%     
    %=====================================================================
    % my version start
    %=====================================================================
    if desperate ~=0
        KbWait(-3);
%         KbPressWait(deviceNumber, 2);
    end
    desperate = 1;
   % KbQueueCreate
     [keyIsDown, secs, keyCode, deltaSecs] = KbCheck;
     
     % important
%      if keyIsDown == 1
%          pause(.2)
%      end
     


%     for frameframe = numFrames*10
%    
%      %[pressed, firstPress, firstRelease, lastPress, lastRelease] = KbQueueCheck(leftKey);
%     end
%     
    % pauses for a fraction of a second when a key is pressed so that I
    % only record one keypress at a time instead of accidentally recording
    % multiple keypresses since it's refreshing 60 times per second (more
    % or less, given the frame refresh rate). Unless they are holding it
    % down, if they hold it down for more than half a second, allow them to
    % scroll quickly
    
    keyscroll_helper = keyscroll_helper + 1;
    keyscroll_tracker(keyscroll_helper) = 0;
    if keyIsDown == 1 && keyCode(five) ==0 && keyCode(fivefive) == 0
        keyscroll_tracker(keyscroll_helper) = 1;
        % if the key has been held for the last half a second, as evidenced
        % by the last 30 keyscroll_tracker columns (if the frame rate is
        % 60) being filled with 1's, then allow the scroller to go much
        % quicker. Only starts after half a second
        if length(keyscroll_tracker) < numFrames/2
            pause(.1)
        else
            % see if the during the last 5 refreshes the key was down, will
            % mean that the last .1*10 = 1 second had the key pressed
            % (due to the pause)
            if sum(keyscroll_tracker(round(end - 10): end)) < length(keyscroll_tracker(round(end - 10): end))/2
                pause(.1)
            end
        end
    end
    
    
    if keyCode(leftKey) == 1 && position > 0
        x = x - (length(scaleRange)/20);
        position = position - 1;
    elseif keyCode(rightKey) == 1 && position < 20
        x = x + (length(scaleRange)/20);
        position = position + 1;
    end
    
    
    SetMouse(round(x), round(rect(4)*scalaPosition), screenPointer, 1);
    
    
    %=====================================================================
    % my version end
    %=====================================================================
    
    % Draw image if provided
    if drawImage == 1
         Screen('DrawTexture', screenPointer, stimuli,[] , rectImage, 0);
    end
    
    % Drawing the question as text
    DrawFormattedText(screenPointer, question, 'center', rect(4)*(scalaPosition - 0.1)); 
    

    % drawing all the 21 points
    offer = 0;
    for tickem = 1:length(all_ticks)
        Screen('DrawLine', screenPointer, scaleColor, all_ticks(tickem, 1), all_ticks(tickem, 2), all_ticks(tickem, 3), all_ticks(tickem, 4), width); % all ticks
        DrawFormattedText(screenPointer, num2str(offer), all_ticks(tickem, 1) - textBounds(1, 3)/2,  rect(4)*scalaPosition+40, [],[],[],[],[],[],[]); % Left point
        offer = offer+1;
    end
    
    Screen('DrawLine', screenPointer, scaleColor, horzLine(1), horzLine(2), horzLine(3), horzLine(4), width);     % Horizontal line
    
    % The slider
    Screen('DrawLine', screenPointer, sliderColor, x, rect(4)*scalaPosition - lineLength*2, x, rect(4)*scalaPosition  + lineLength*2, width*2);
    
    % Caculates position
%     position          = round((x)-mean(scaleRange));           % Shift the x value according to the new scale
%     position          = (position/max(scaleRangeShifted))*100; % Converts the value to percentage
    
    % Display position
    if displayPos
        DrawFormattedText(screenPointer, num2str(round(position)), 'center', rect(4)*(scalaPosition + 0.05)); 
    end
    
    % Flip screen
    onsetStimulus = Screen('Flip', screenPointer);
  
     
    % Check if answer has been given
    if strcmp(device, 'mouse')
        secs = GetSecs;
        if buttons(responseKey) == 1
            answer = 1;
        end
    elseif strcmp(device, 'keyboard')
        [keyIsDown, secs, keyCode] = KbCheck;
        if keyCode(responseKey) == 1
            answer = 1;
            RT = secs - t0;
        end
    else
        error('Unknown device');
    end
    
    
    
    
    % Abort if answer takes too long
    if secs - t0 > aborttime 
        RT = secs - t0;
        break
    end
end

% freeze on the choice for a half second plus a half second jitter
% % still_freeze = 0;
% % % while still_freeze == 0
for frame = 1:numFrames/2 + (numFrames/2*rand)
    % drawing the tick where the chosen offer is
    
    for tickem = 1:length(all_ticks)
        if tickem == position+1
            Screen('TextStyle', screenPointer, 1);
            
            Screen('DrawLine', screenPointer, scaleColor, all_ticks(tickem, 1), all_ticks(tickem, 2), all_ticks(tickem, 3), all_ticks(tickem, 4), width); % all ticks
            DrawFormattedText(screenPointer, num2str(position), all_ticks(tickem, 1) - textBounds(1, 3)/2,  rect(4)*scalaPosition+40, [],[],[],[],[],[],[]);
            
            Screen('TextStyle', screenPointer, 0);
        elseif tickem == 1
            Screen('DrawLine', screenPointer, scaleColor, all_ticks(tickem, 1), all_ticks(tickem, 2), all_ticks(tickem, 3), all_ticks(tickem, 4), width); % all ticks
        elseif tickem == length(all_ticks)
            Screen('DrawLine', screenPointer, scaleColor, all_ticks(tickem, 1), all_ticks(tickem, 2), all_ticks(tickem, 3), all_ticks(tickem, 4), width); % all ticks
        end
        
    end
    
    % Display position
    
% %     DrawFormattedText(screenPointer, num2str(round(position)), 'center', rect(4)*(scalaPosition + 0.05)); 
    
    Screen('DrawLine', screenPointer, scaleColor, horzLine(1), horzLine(2), horzLine(3), horzLine(4), width);     % Horizontal line
    
    % The slider
    Screen('DrawLine', screenPointer, sliderColor, x, rect(4)*scalaPosition - lineLength*2, x, rect(4)*scalaPosition  + lineLength*2, width*2);
    onsetStimulus = Screen('Flip', screenPointer);
    
% %     pause(1)
    
% %     still_freeze = 1;
end



%% Wating that all keys are released and delete texture
KbReleaseWait; %Keyboard
% KbReleaseWait(1); %Mouse
if drawImage == 1
    Screen('Close', stimuli);
end
%% Calculating the rection time and the position
% this is actually calculated earlier in the script now
% RT                = (secs - t0)*1000;                                          % converting RT to millisecond
% % RT                = (secs - t0);                                          % converting RT to millisecond
end