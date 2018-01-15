function [done] = do_task_perceptconf(task,nsub,nsession,DT)


% Initialize

if task == 1
    resultname=strcat('ConfConIncentive_Sub',num2str(nsub),'_Session',num2str(nsession));
    calibname = strcat('PerceptCalibCon_Sub',num2str(nsub),'_Session1');
elseif task == 2
    resultname=strcat('ConfConShort_Sub',num2str(nsub),'_Session',num2str(nsession));
    calibname = strcat('PerceptCalibCon_Sub',num2str(nsub),'_Session2');
elseif task == 3
    resultname=strcat('ConfConHybrid_Sub',num2str(nsub),'_Session',num2str(nsession));
    calibname = strcat('PerceptCalibCon_Sub',num2str(nsub),'_Session2');
end


calib = load(calibname);


% load incentive pictures
%=========================%

%loadsound('big_gain.wav',1); loadsound('small_gain.wav',2); loadsound('big_loss.wav',3); loadsound('small_loss.wav',4); loadsound('neutral.wav',5);

cgloadbmp(2,'euroGain100.bmp'); cgloadbmp(3,'centGain100.bmp'); cgloadbmp(4,'euroLoss100.bmp'); cgloadbmp(5,'centLoss100.bmp');
cgloadbmp(6,'NoeuroGain100.bmp'); cgloadbmp(7,'NocentGain100.bmp'); cgloadbmp(8,'NoeuroLoss100.bmp'); cgloadbmp(9,'NocentLoss100.bmp');

Inc = {'You can win 1 euro!';'You can win 10 cents!';'You can lose 1 euro!';'You can lose 10 cents!'};
FBwin = {'You win 1 euro!';'You win 10 cents!';'You don''t lose 1 euro!';'You don''t lose 10 cents!'};
FBloss = {'You don''t win 1 euro!';'You don''t win 10 cents!';'You lose 1 euro!';'You lose 10 cents!'};

inc_mat = [2,3,4,5;6,7,4,5;2,3,8,9]';
sound_mat = [5,5,3,4;1,2,5,5]';

% parametrize gabor stuff
%==========================%
G.pos = 200;                               % image position
G.imSize = 200;                            % image size: n X n
G.lamda = 20.*(G.imSize/100);              % wavelength (number of pixels per cycle)
G.sigma = 15.*(G.imSize/100);              % gaussian standard deviation in pixels, if 0 then don't apply gaussian
G.phase = .25;                             % phase (0 -> 1)
G.trim = .001;                             % trim off edges of gaussian with values < trim
G.theta = 15;                              % grating orientation (in degrees)
G.con = .1;                                 % diminish the contrast
G1 = G; G2 = G;             % attribute properties to G1 and G2
G1.pos = -G.pos;            % change side of G1: left gabor


% Parametrize lottery structures
%==============================%
L.coord = [0,-100];                         % lottery position
L.width = 250;                              % lottery width
L.ptext_coord = [0,-100 + L.width/2+30];    % lottery probabilitry position
L.winning_color = [0,1,0];                  % color of the winning area (RGB)
L.losing_color = [1,0,0];                   % color of the losing area (RGB)
L.spin = [];                                % spinning
L.win = [];                                 % winning
L.FBwin = FBwin;                            % Feedback for winning
L.FBloss = FBloss;                          % Feedback for losing
L.inc = [];                                 % Incentive category

% Parametrize clock structures
%==============================%
C.coord = [0,-100];                         % lottery position
C.width = 250;                              % lottery width
C.text_coord = [0,-100 + L.width/2+30];    % lottery txt position
C.color = .8*[1,1,1];                  % color of the winning area (RGB)
C.spin = [];                                % spinning
C.win=[];                                   % winning
C.FBwin = FBwin;                            % Feedback for winning
C.FBloss = FBloss;                          % Feedback for losing
C.inc = [];                                 % Incentive category


% parametrize confidence rating scale
%======================================%
K.coord = [0,-100];                            % scale position
K.width = 800;                              % scale width (in pixels)
K.color = [0,0,0];                          % scale color (RGB)
K.minmax = [50,100];                        % scale numbers mins and max
K.Mgrad = 6;                                % scale Major Ticks
K.mgrad = 11;                               % scale minor Ticks
K.curs = 0;                                 % cursor position
K.curscol = [1,1,0];                        % cursor color

% Get design features
%=====================%
fixation_time = 750;
display_time = DT;
response_time = 100;
incentive_time = 1000;
comp1_time = 1200;
comp2_time = 1200;
lot_time = 750;
result_time = 1000;

% Get Keys
%===========%
leftkey = 97;
rightkey = 98;
spacebar = 71;
escapekey = 52;

% Get Key variables, according to the calib
%============================================%
pc = repmat((.15:.05:.85)',2,1);
DV = 1./calib.GLM.b(2).*(log(pc./(1-pc))-calib.GLM.b(1));% here DV is G2con - G1con
mDV = repmat(DV,4,1); mpc = repmat(pc,4,1);
bb = ones(length(pc),1); ninc = [bb;2*bb;3*bb;4*bb];

m = NaN(length(mDV),1);
m(mDV>0)=0;m(mDV<0)=-mDV(mDV<0);

G1val = m + .5*rand(length(mDV),1);
G2val = G1val + mDV;

% Randomize key variables
%==========================%
kk = randperm(4*length(pc));
G1val = G1val(kk); G2val = G2val(kk);
mDV = mDV(kk); mpc = mpc(kk);
ninc = ninc(kk);

% Compute decision variables to check the design
%================================================%
check = [G1val,G2val,mDV,mpc];
XX = 1./(1+exp(-(calib.GLM.b(1)+calib.GLM.b(2).*(check(:,2)-check(:,1)))));
check = [check,XX,ninc];

% Get Data to Save
%=================%
n_trial = size(G1val,1);

trials = (1:n_trial)';
time0 = NaN(n_trial,1);
time1 = NaN(n_trial,1);
time2 = NaN(n_trial,1);
time3 = NaN(n_trial,1);
time4 = NaN(n_trial,1);
resp = NaN(n_trial,1);
time5 = NaN(n_trial,1);
time6 = NaN(n_trial,1);
initpos = floor(5*rand(n_trial,1)-2);
conf = NaN(n_trial,1);
lot = (50 + floor(51*rand(n_trial,1)))./100;
time7 = NaN(n_trial,1);
time8 = NaN(n_trial,1);
rand_trial_p = floor(rand(n_trial,1)*100+1)/100;
lot_win = NaN(n_trial,1);
true_ans = NaN(n_trial,1);
resp_win = NaN(n_trial,1);


% Ready to start
%=================%
setforecolour(0,0,0);
settextstyle('Arial', 32 )

preparestring('Your goal is to indicate, at each trial,',1,0,100);
preparestring('which of the 2 gabor patches',1,0,50);
preparestring('has the highest constrast,',1,0,0);
preparestring('and to indicate your confidence in your answer.',1,0,-50);
preparestring('Press the spacebar to continue',1,0,-200);
drawpict(1);
[keydown,~,~]=waitkeydown(inf,[spacebar,escapekey]);
if keydown == escapekey
    stop_cogent
end
clearpict(1);

preparestring('Ready ?',1,0,50);
preparestring('Press the spacebar to start',1,0,-50);
drawpict(1);
[keydown,~,~]=waitkeydown(inf,[spacebar,escapekey]);
if keydown == escapekey
    stop_cogent
end
clearpict(1);


% Trial Loop
%============%
T0 = time;
for k_trial = 1:n_trial
    
    % Get gabors parameters for this trial
    G1.theta = ceil(360*rand());
    G2.theta = ceil(360*rand());
    
    G1.con = G1val(k_trial);
    G2.con = G2val(k_trial);
    
    [~,im1] = gaborFn_ML(G1.imSize, G1.lamda, G1.sigma, G1.theta, G1.phase, G1.trim, G1.con);  % gabor function
    [~,im2] = gaborFn_ML(G2.imSize, G2.lamda, G2.sigma, G2.theta, G2.phase, G2.trim, G2.con);  % gabor function
    
    % Fixation
    setforecolour(0,0,0);
    settextstyle('Arial', 55 )
    preparestring('+',1,0,0);
    time0(k_trial) = drawpict(1) - T0;
    waituntil(T0 + time0(k_trial) + fixation_time);
    
    % Display gabors
    preparepict(im1,1,G1.pos,0);
    preparepict(im2,1,G2.pos,0);
    time1(k_trial) = drawpict(1) - T0;
    waituntil(T0 + time1(k_trial) + display_time);
    clearpict(1,.5,.5,.5);
    
    % response self paced
    clearkeys;
    settextstyle('Arial', 55 )
    preparestring('?',1,0,0);
    preparestring('<',1,-30,0);
    preparestring('>',1,30,0);
    time2(k_trial) = drawpict(1) - T0;
    [keydown,timedown,numberdown] = waitkeydown(Inf,[leftkey,rightkey,escapekey]);
    time3(k_trial) = timedown - T0;
    resp(k_trial) = (keydown == leftkey)*2-1;
    
    % response selection
    setforecolour(1,0,0)
    settextstyle('Arial', 55 )
    if keydown == leftkey
        preparestring('<',1,-30,0);
    elseif keydown == rightkey
        preparestring('>',1,30,0);
    elseif keydown == escapekey
        stop_cogent
    end
    drawpict(1);
    waituntil(T0 + time3(k_trial) + response_time);
    clearpict(1,.5,.5,.5);
    time4(k_trial) = drawpict(1) - T0;
    
    
    % Incentive
    cgflip(.5,.5,.5)
    cgfont('Arial',38);cgpencol(0,0,0);
    cgtext(Inc{ninc(k_trial)},0,0)
    cgdrawsprite(inc_mat(ninc(k_trial),1),0,200)
    cgflip(.5,.5,.5)
    waituntil(T0 + time4(k_trial) + incentive_time);
    
    % confidence rating
    clearkeys
    readkeys
    step = round(K.width./(K.mgrad-1));
    % initialise params
    K.lot = [];
    K.type =[];
    K.curs = initpos(k_trial).*step;
    % reday to rate
    stop_signal = 0;
    time5(k_trial) = time - T0;
    
    while stop_signal == 0;
        viz_gauge(K)
        cgdrawsprite(inc_mat(ninc(k_trial),1),0,200)
        cgflip(.5,.5,.5)
        [keydown,timedown,~] = waitkeydown(inf, [leftkey,rightkey,spacebar,escapekey]);
        if keydown == spacebar
            stop_signal = 1;
        elseif keydown == leftkey
            x = max([K.coord(1)-K.width./2,K.curs-step]);
            K.curs = x(1);
        elseif keydown == rightkey
            x = min([K.coord(1)+K.width./2,K.curs+step]);
            K.curs = x(1);
        elseif keydown == escapekey
            stop_cogent
        end
    end
    conf(k_trial) = (50+(5+K.curs./step)*5)./100;
    true_ans(k_trial) = (G1.con>G2.con)*2-1;
    resp_win(k_trial) = resp(k_trial) == true_ans(k_trial);
    
    % confidence vs lottery determines gain
    time6(k_trial) = time - T0;
    if task == 1 || task ==3
        K.conf = conf(k_trial); K.lot = lot(k_trial);
        cgdrawsprite(inc_mat(ninc(k_trial),1),0,200)
        viz_gauge(K)
        cgflip(.5,.5,.5)
        
        if task == 1
            waituntil(T0 + time6(k_trial) + comp1_time)
        end
        
        K.type = 2*(K.conf>=K.lot)-1;
        cgdrawsprite(inc_mat(ninc(k_trial),1),0,200)
        viz_gauge(K)
        cgflip(.5,.5,.5)
        waituntil(T0 + time6(k_trial) + comp1_time + (task==1).*comp2_time)
    end
    
    %----------------%
    time7(k_trial) = time - T0;
    % Case 1 : answer determines gain
    if conf(k_trial) >= lot(k_trial)
        
        C.win = []; C.spin = []; C.inc = ninc(k_trial);
        if task ==1
            % determine gain
            % Display lottery
            viz_clock(C)
            cgdrawsprite(inc_mat(ninc(k_trial),1),0,200)
            cgflip(.5,.5,.5)
            waituntil(time7(k_trial)+ T0 + lot_time)
            % Spin lottery
            tclock = 750 + round(500*rand());
            ntick = 60;
            tunit = tclock/(ntick+1);
            instant_wheel = 0;                                          % first position of the spining wheel
            for k_ang = 0:ntick;
                t2 = time;
                clock_spin = k_ang*2*pi/ntick;
                cgdrawsprite(inc_mat(ninc(k_trial),1),0,200)
                C.spin = clock_spin;
                viz_clock(C)
                cgflip(.5,.5,.5);
                waituntil(t2 + tunit)                                        % wait 25 ms after the spin
            end
        end
        
        if task == 1 || task ==3
            % Display result
            time8(k_trial) = time - T0;
            C.win = 3+resp_win(k_trial);
            cgdrawsprite(inc_mat(ninc(k_trial),resp_win(k_trial)+2),0,200)
            viz_outcome(C)
            cgflip(.5,.5,.5)
            %playsound(sound_mat(ninc(k_trial),resp_win(k_trial)+1));
            waituntil(T0 + time8(k_trial) + result_time)                                    % display result for 1000 ms after spin
        end
        
        % Case 2 : lottery determines gain
    elseif conf(k_trial) < lot(k_trial)
        L.proba = lot(k_trial);
        L.win = []; L.spin = []; L.inc = ninc(k_trial);
        
        if task ==1
            % determine gain
            % Display lottery
            cgdrawsprite(inc_mat(ninc(k_trial),1),0,200)
            viz_wheels(L)
            cgflip(.5,.5,.5)
            waituntil(time7(k_trial)+ T0 + lot_time)
            % Spin lottery
            ang = 1 + rand_trial_p(k_trial);
            instant_wheel = 0;                                          % first position of the spining wheel
            while abs(instant_wheel - ang) > 0.002
                t4 = time;
                lottery_spin = instant_wheel * 2 * pi;
                cgdrawsprite(inc_mat(ninc(k_trial),1),0,200)
                viz_wheels(L)
                L.spin = lottery_spin;
                cgflip(.5,.5,.5)
                if abs(instant_wheel - ang) < 0.75                       % Determine angluar speed and slow down sppinning when close to target
                    instant_angle = instant_angle * 0.96;
                else instant_angle = 0.035;
                end
                instant_wheel = instant_wheel + (instant_wheel < ang) * instant_angle - (ang < instant_wheel) * instant_angle;
                waituntil(t4 + 25)                                        % wait 25 ms after the spin
            end
        end
        
        if task == 1 || task ==3
            % Display result
            time8(k_trial) = time - T0;
            lot_win(k_trial) = double(rand_trial_p(k_trial) < L.proba);
            L.win = 1+lot_win(k_trial);
            cgdrawsprite(inc_mat(ninc(k_trial),lot_win(k_trial)+2),0,200)
            viz_outcome(L)
            cgflip(.5,.5,.5)
            %playsound(sound_mat(ninc(k_trial),lot_win(k_trial)+1));
            waituntil(T0 + time8(k_trial) + result_time)
        end% display result for 1000 ms after spin
    end
    
    
    clearpict(1,.5,.5,.5)
    
    if task ==2
        if mod(k_trial,20)==0;
            setforecolour(0,0,0);
            settextstyle('Arial', 48 )
            cperf = mean(resp_win(k_trial-19:k_trial));
            preparestring('Your performance in the past 20 trials is',1,0,50);
            preparestring(strcat(num2str(cperf*100),'% accuracy'),1,0,0);
            preparestring('Press spacebar to continue',1,0,-100);
            drawpict(1);
            waitkeydown(inf,spacebar);
        end
        clearpict(1)
        drawpict(1)
        
    end
end

datatime.all = [time0,time1,time2,time3,time4,time5,time6,time7,time8];
datatime.t0 = T0;
data = [trials,G1val,G2val,resp,ninc,initpos,conf,lot,rand_trial_p,lot_win,true_ans,resp_win];

save(resultname,'data','datatime');

clearpict(1,.5,.5,.5); drawpict(1);
setforecolour(0,0,0);
settextstyle('Arial', 32 )
preparestring('This task is now finished',1,0,50);
preparestring('You have earned X euros extra!,1,0,0');
preparestring('Please press the spacebar to continue.',1,0,-50);
drawpict(1);
[keydown,~,~]=waitkeydown(inf,spacebar);
if keydown == spacebar
    clearpict(1);cgflip;drawpict(1); clc;
    done =1;
end

