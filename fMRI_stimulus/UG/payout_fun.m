% function [sub_payout] = payout_fun
% this function calculates the final extra payout of the participant

% For the probability matching task, I am going to make the trail worth
% 2.50 (or 20 MU), since that's the maximum that they can make in the UG.
% If I understand it correctly, the probability matching is an all or
% nothing gamble, in which subjects increase their chances of getting the
% prize the closer they are to the true probability

% Script must be called the directory where the script is

tmp_sub_num = input('Enter the subject number: \n');

% row 1 is the payout from the UG
% row 2 is the payout from the probability matching task
% row 3 is the payout from the dictator game
% row 4 is the payout from the Eckel-Grossman risk lottery
% importantly, all these values are in euro, this function converts from MU
% to euro
sub_payout_mat = NaN(1, 4);

%% Getting all the data indices
if tmp_sub_num > 99
    sub_num = ['sub' num2str(tmp_sub_num)];
elseif tmp_sub_num > 9 && tmp_sub_num < 100
    sub_num = ['sub0' num2str(tmp_sub_num)];
elseif tmp_sub_num < 10
    sub_num = ['sub00' num2str(tmp_sub_num)];
end

sub_files     = dir([pwd, filesep, 'data', filesep,  '*', sub_num, '*']);

dict_risk_idx = NaN(length(sub_files), 1);
prob_idx      = NaN(length(sub_files), 1);
UG_idx        = NaN(length(sub_files), 1);

for ii = 1:length(sub_files)
    if contains(sub_files(ii).name, 'dict_risk') == 1
        dict_risk_idx(ii) = ii;
    elseif contains(sub_files(ii).name, 'prob') == 1
        prob_idx(ii) = ii;
    elseif contains(sub_files(ii).name, 'soc') == 1
        UG_idx(ii) = ii;
    end
end

dict_risk_idx(any(isnan(dict_risk_idx), 2), :) = [];
prob_idx(any(isnan(prob_idx), 2), :) = [];
UG_idx(any(isnan(UG_idx), 2), :) = [];

% if there's two dict_risk_idx it means that one has a time and date stamp,
% and that one came aftet the first so that's the one that I'll use
if length(dict_risk_idx) > 1
    dict_risk_idx = dict_risk_idx(2);
end
% same is true of prob_idx
if length(prob_idx) > 1
    prob_idx = prob_idx(2);
end

% if there's more than one ultimatum game file with time-stamps, I'm just
% going to leave them in the mix, the trial will still be chosen randomly 
UG_idx = Shuffle(UG_idx); % Shuffle up the indices 
UG_idx = UG_idx(1);       % picking the first one after shuffling is picking at random

%% UG
% Selecting the UG winning trial and convering from MU to euro
load(['data', filesep, sub_files(UG_idx).name]);           % loads in the UG .mat file, which includes sub_data
all_rows = Shuffle(1:length(sub_data));

sub_payout_mat(1) = sub_data(all_rows(1), 3) * .125; % row 3 is the payoff row, .125 is the MU conversion rate


%% Probability matching
% selecting the probability matching trial, determing if it's a victory or
% not (if it's a victory subject receives 2.50)
load(['data', filesep, sub_files(prob_idx).name]);
all_rows = Shuffle(1:length(sub_data));

lotto = datasample(0:100, 1); % sample the lottery from the same selection the subjects had

prob_row = sub_data(all_rows(1), [1,2,4]);

if prob_row(3) == 0
    % endowment of 0 true probabilities
    true_probs =  [0.0707, 0.1032, 0.1482, 0.2083, 0.2846, 0.3757, 0.4764, 0.5791, 0.6754, 0.7589, 0.8264, 0.8780, 0.9159, 0.9427, 0.9614, 0.9741, 0.9827, 0.9885, 0.9924, 0.9949, 0.9967];
elseif prob_row(3) == 1
    % endowment of 10 true probabilities
    true_probs = [0.2844, 0.3839, 0.4942, 0.6051, 0.7061, 0.7902, 0.8552, 0.9026, 0.9356, 0.9579, 0.9728, 0.9825, 0.9887, 0.9928, 0.9954, 0.9971, 0.9981, 0.9988, 0.9992, 0.9995, 0.9997];
elseif prob_row(3) == 2
    true_probs = [0.8268, 0.8855, 0.9261, 0.9530, 0.9705, 0.9816, 0.9885, 0.9929, 0.9956, 0.9973, 0.9983, 0.9990, 0.9994, 0.9996, 0.9998, 0.9998, 0.9999, 0.9999, 1.0000, 1.0000, 1.0000];
end

if prob_row(2) > lotto
%     win = rand(1) < abs((prob_row(2) * .01) - true_probs(prob_row(1)));
    win = rand(1) < true_probs(prob_row(1));
else
    win = rand(1) < (lotto * .01);
end

if win
   sub_payout_mat(2) = 2.5; % if they win, they get 2.50 euro
else 
   sub_payout_mat(2) = 0; % if they lose, they get 0 euro
end


%% Dictator
load(['data', filesep, sub_files(dict_risk_idx).name]);           % loads in the Dict .mat file, which includes sub_data
all_rows = Shuffle(1:length(sub_data));

sub_payout_mat(3) = (20 - sub_data(all_rows(1), 3)) * .125; % row 1 is the offer, .125 the MU conversion


%% Risk
% this is contained in the same .mat file as the dicator
% 0 = low, 1 = high, chosen at %50 probability
highlow = datasample(0:1, 1);
if risk_gamble == 1
    sub_payout_mat(4) = 28*.125;
elseif risk_gamble == 2
    if highlow == 0
        sub_payout_mat(4) = 24*.125;
    elseif highlow == 1
        sub_payout_mat(4) = 36*.125;
    end
elseif risk_gamble == 3
    if highlow == 0
        sub_payout_mat(4) = 20*.125;
    elseif highlow == 1
        sub_payout_mat(4) = 44*.125;
    end
elseif risk_gamble == 4
    if highlow == 0
        sub_payout_mat(4) = 16*.125;
    elseif highlow == 1
        sub_payout_mat(4) = 52*.125;
    end
elseif risk_gamble == 5
    if highlow == 0
        sub_payout_mat(4) = 12*.125;
    elseif highlow == 1
        sub_payout_mat(4) = 60*.125;
    end
elseif risk_gamble == 6
    if highlow == 0
        sub_payout_mat(4) = 2*.125;
    elseif highlow == 1
        sub_payout_mat(4) = 70*.125;
    end
end

%% Final
sub_payout = sum(sub_payout_mat);

fprintf('%f\n', sub_payout);
