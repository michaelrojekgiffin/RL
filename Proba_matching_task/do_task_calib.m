
function [done] = do_task_calib(nsub,nsess,DT)

resultname=strcat('PerceptCalibCon_Sub',num2str(nsub),'_Session',num2str(nsess));



% parametrize gabor stuff
%==========================%
G.pos = 200;                               % image position
G.imSize = 200;                            % image size: n X n
G.lamda = 20.*(G.imSize/100);              % wavelength (number of pixels per cycle)
G.sigma = 15.*(G.imSize/100);              % gaussian standard deviation in pixels, if 0 then don't apply gaussian
G.phase = .25;                             % phase (0 -> 1)
G.trim = .001;                             % trim off edges of gaussian with values < trim
G.theta = 15;                              % grating orientation (in degrees)
G.con = .05;                                 % diminish the contrast
G1 = G; G2 = G;             % attribute properties to G1 and G2
G1.pos = -G.pos;            % change side of G1: left gabor

% Get design features
%=====================%
n_block = 11;
n_trialpblock = 12;
n_trial = n_block.*n_trialpblock;
fixation_time = 750;
display_time = DT;
response_time = 100;

% Get Keys
%===========%
leftkey = 97;
rightkey = 98;
spacebar = 71;
escapekey = 52;

% Get Data to Save
%=================%
trials = (1:n_trial)';
G1val = NaN(n_trial,1);
G2val = NaN(n_trial,1);
time0 = NaN(n_trial,1);
time1 = NaN(n_trial,1);
time2 = NaN(n_trial,1);
time3 = NaN(n_trial,1);
time4 = NaN(n_trial,1);
resp = NaN(n_trial,1);
true_ans = NaN(n_trial,1);
resp_win = NaN(n_trial,1);
time5 = NaN(n_trial,1);
time6 = NaN(n_trial,1);


block_perf = NaN(n_block,1);
mult =  NaN(n_block,1);

% Ready to start
%=================%
setforecolour(0,0,0);
settextstyle('Arial', 32 )

preparestring('Your goal is to indicate, at each trial,',1,0,100);
preparestring('which of the 2 gabor patches',1,0,50);
preparestring('has the highest contrast.',1,0,0);
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
mult(1) = 1;
k_trial = 0;
for k_block = 1:n_block
    
    k_in = (k_block-1)*n_trialpblock + 1;
    k_out = k_block*n_trialpblock;
    
    
    kk = randperm(n_trialpblock);
    zz = linspace(0,.5,n_trialpblock/2);
    xx = .5*rand(1,n_trialpblock/2);
    g1temp = [xx,xx+mult(k_block).*zz];
    g2temp = [xx+mult(k_block).*zz,xx];
    
    G1val(k_in:k_out,1) = g1temp(kk);
    G2val(k_in:k_out,1) = g2temp(kk);
    
    
    block_perf(k_block) = 0;
    
    for k = 1:n_trialpblock
        
        k_trial = k_trial +1;
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
        time2(k_trial) = drawpict(1) - T0;
        
        % response self paced
        clearkeys;
        settextstyle('Arial', 55 )
        preparestring('?',1,0,0);
        preparestring('<',1,-30,0);
        preparestring('>',1,30,0);
        time3(k_trial) = drawpict(1) - T0;
        [keydown,timedown,numberdown] = waitkeydown(Inf,[leftkey,rightkey,escapekey]);
        time4(k_trial) = timedown - T0;
        resp(k_trial) = (keydown == leftkey)*2-1;
        
        true_ans(k_trial) = (G1.con>G2.con)*2-1;
        resp_win(k_trial) = resp(k_trial) == true_ans(k_trial);
        
        block_perf(k_block) = block_perf(k_block) + resp_win(k_trial)./n_trialpblock;
        
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
        time5(k_trial) = drawpict(1) - T0;
        waituntil(T0 + time5(k_trial) + response_time);
        clearpict(1,.5,.5,.5);
        time6(k_trial) = drawpict(1) - T0;
        cgflip(.5,.5,.5)
        
    end
    
    if k_block == 1
        mult(k_block+1) = mult(k_block);
    elseif k_block < n_block && k_block>1
        if block_perf(k_block) < .7
            mult(k_block+1) = min([1,mult(k_block).*1.2]);
        elseif block_perf(k_block) > .7
            mult(k_block+1) = max([.1,mult(k_block).*.8]);
        elseif block_perf(k_block) == .7
            mult(k_block+1) = mult(k_block);
        end
    end
    
end

datablock = [mult,block_perf];
datatime.all = [time0,time1,time2,time3,time4,time5,time6];
datatime.t0 = T0;
datafull = [trials,G1val,G2val,resp,true_ans,resp_win];
data = datafull(n_trialpblock+1:end,:);
[GLM.b,~,GLM.stats]=glmfit(data(:,3)-data(:,2),(data(:,4)+1)./2,'binomial');

save(resultname,'data','datatime','datafull','datablock','GLM');

clearpict(1,.5,.5,.5); drawpict(1);
setforecolour(0,0,0);
settextstyle('Arial', 32 )
preparestring('This task is now finished',1,0,50);
preparestring('Please press spacebar to continue with the main task.',1,0,-50);
drawpict(1);
[keydown,~,~]=waitkeydown(inf,spacebar);
if keydown == spacebar
    clearpict(1);cgflip;drawpict(1); clc;
    done =1;
end
