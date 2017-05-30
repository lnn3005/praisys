%%% Simulation for RCPSP

close all
clear
clc

%% Get input

filename1 = 'input-task-43.txt';
filename2 = 'input-relation4.txt';
filename3 = 'input-task-43-result.txt';

task = dlmread(filename1); 
rel = dlmread(filename2);

% Vector of Initial state 

% Vector of Optimal schedule