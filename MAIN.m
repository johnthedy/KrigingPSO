clc;clear;tic;close all
format long
%% Input Explanation
% probtype : Choose problem type from 1 to 14, if user have own problem, then edit it on G.m and problem.m
% uxdoe and lxdoe : upper and lower bound of standard normal space that is considered.
% beta1 and beta2 : initial outter and inner radius in standard normal space, it is suggested to use 8 and 3 respectively.
% Kriginginitial : initial sample size to build Kriging, just use as small as possible, this study suggest to use 5 x number of random variables depends on problem complexity, 
% PfCOV : Monte Carlo Simulation Coefficient of Variation, expected COV on final failure probability evaluation.
% Ulimit and REIFlimit : At every iteration, best particle is obtained from PSO, its Ufunction and REIFfunction value is evaluated, limit value for Ufunction and REIFfunction value is set so that next level is given if this threshold value get crossed.
% plotfig : put 1 if you want to plot figure, only applicable for 2D problem.

% This code require dace Kriging algorithm to run.
% Dace Kriging algorithm is solely taken from 
% Nielsen, H. B. (Author), Lophaven, S. N. (Author), & Søndergaard, J. (Author). (2002). 
% DACE - A Matlab Kriging Toolbox. Computer programme, Informatics and Mathematical Modelling, Technical University of Denmark, DTU. http://www2.imm.dtu.dk/pubdb/p.php?1460
%
% Many thanks to Nielsen, H. B. (Author), Lophaven, S. N. (Author), & Søndergaard, J. (Author). (2002). for providing DACE Kriging algorithm
%% Kriging + PSO + HHS

probtype=4;
uxdoe=8;
lxdoe=-8;
beta1=8;
beta2=3;
Kriginginitial=10;
PfCOV=0.02;
Ulimit=1.4;
REIFlimit=0.05;
plotfig=1;

currentpath=sprintf('%s\\%s',pwd,'dace');
addpath(currentpath)
[Pf,FE]=KrigingPSO(probtype,uxdoe,lxdoe,beta1,beta2,Kriginginitial,PfCOV,Ulimit,REIFlimit,plotfig);
