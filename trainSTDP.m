%% Train network with given speech patterns for STDP rules
% assume following steps done
% -> audio-> Lyon Auditory filter model (or some model)->BSA (or other)
% ->1KHz sampling -> Network input in *.mat file
clear all; close all; clc;
set(0,'DefaultFigureWindowStyle','docked');
%% Get liquid network lattice
