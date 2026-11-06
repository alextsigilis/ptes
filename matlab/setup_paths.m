% setup_paths.m
% Add every folder of the MATLAB code base to the search path.
% Run once per session from anywhere:  >> run('matlab/setup_paths.m')
root = fileparts(mfilename('fullpath'));
addpath(genpath(root));
fprintf("Added %s (and subfolders) to the MATLAB path.\n", root);
