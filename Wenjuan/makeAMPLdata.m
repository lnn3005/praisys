
close all
clear
clc

%% Read input

filename1 = 'input-task-43.txt';
filename2 = 'input-relation4.txt';

task = dlmread(filename1); 
rel = dlmread(filename2);

%% Amount of resources available per period - Hard-coded
if filename1 == 'input-task-42.txt' 
    availResources = [8 6 6];
end

if filename1 == 'input-task-43.txt'
    availResources = [6];
end

%% Parameters

numRows = size(task,1);
numCols = size(task,2);

% The number of tasks is the number of rows
numTasks = numRows;

% The number of modes is 2 
% Hard-coded
numModes = 2;

% The number of resource type 
numResourcesTypes = (numCols - 3)/2;

% Resources matrix
resources = zeros(numTasks,numModes,numResourcesTypes);
start_col = 1;
end_col = 1;
for i=1:numModes
    % Start column of mode i  
    start_col = end_col + 2;
    % End column of mode i
    end_col = start_col + numResourcesTypes - 1;
    resources(:,i,:) = task(:,start_col:end_col);
end

% Duration is represented column 2 for mode 1 and
% columnn 3+numResourcesTypes for mode 2
duration = task(:,[2 3+numResourcesTypes]);

% Set horizon length
horizonLength = sum(max(duration,[],2));

% The set of tasks is 1 2 ... numTasks
taskSet = 1:numTasks;

% Set earliest start date of all tasks to be 1
earliest = ones(1,numTasks);

% Set latest start date of all tasks to be the last period
latest = zeros(1,numTasks) + horizonLength;


%% Start AMPL data file 

% Create a subfolder called 'AMPL data files'. A new data file is created
% inside this folder each time the code is run.

dir = 'AMPL data files';
filename = strcat(datestr(now,'ddmmmyy_HHMM'),'.dat');

if ~exist(dir,'dir')  
    mkdir(dir);
end


%% Write to .dat file

% Open data file to be written
fileID = fopen(fullfile(dir,filename),'w');

fprintf(fileID,'data; \n\n');

% Print parameter
writeAMPL(fileID,2,'numTasks',numTasks);
writeAMPL(fileID,2,'horizonLength',horizonLength);
writeAMPL(fileID,2,'numResourcesTypes',numResourcesTypes);
writeAMPL(fileID,2,'Ra',availResources);

% Print precedence set
writeAMPL(fileID,11,'P',rel);

% Print mode sets
for i=1:numTasks
    nameString = ['M[',num2str(i),']'];
    writeAMPL(fileID,1,nameString,1:numModes);
end

% Print earliest start date
writeAMPL(fileID,2,'e',earliest);

% Print latest start date
writeAMPL(fileID,2,'l',latest);

% Print duration
writeAMPL(fileID,2,'d',duration);

% Print resource required
writeAMPL(fileID,23,'w',resources);

fclose(fileID);
