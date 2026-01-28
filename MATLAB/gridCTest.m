
% Script to test gridGenerator
clc; clear; close all;

% Add Grids folder to path
addpath(fullfile(pwd,'Grids'));


% --- Parameters ---
Ns = 41;       % # columns (s-direction)
Nt = 41;       % # rows    (t-direction)
N_basis = 0;  % # basis functions per direction

% --- Grid function handle ---
grid_fun = @makeGridC;  % must return [X,Y] of size Nt x Ns

% --- Run the optimization ---
gridGenerator(grid_fun, Ns, Nt, N_basis);
