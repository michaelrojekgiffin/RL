clear all
close all
clc

% Get path & sub num
%====================%
cur_dir = pwd;
cogent_path = 'D:\Confidence_exp\Cogent2000v1.33\Task\Task.v5';

nsub = input('Subject?');
dt = 150;



% parametrize cogent parameters
%===============================%
addpath(genpath(cogent_path))
config_display(0,3,[0.5,0.5,0.5],[1,1,1], 'Arial', 55, 3);    %easier to work on - window mode
%config_display(1,3,[0 0 0],[1 1 1], 'Arial', 55, 1);
config_keyboard(100,5,'nonexclusive');
config_sound;
start_cogent;


% run tasks
%=============%

% [do_calib1] = do_task_calib(nsub,1,dt);
[do_task1] = do_task_perceptconf(1,nsub,1,dt);
[do_calib2] = do_task_calib(nsub,2,dt);
[do_task2] = do_task_perceptconf(3,nsub,1,dt);
addpath(genpath(vba_path))
rmpath(genpath(vba_path))

% stop
%=============%
stop_cogent
rmpath(genpath(cogent_path))