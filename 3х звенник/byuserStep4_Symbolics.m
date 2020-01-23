close all; clear; %clear classes;
clc; 

%Create user interfase object for SRD
SRD = SRDuserinterface();

%Generate dynamics eq.
ToLinearize = true; 
ToSimplify = true; 
SRD.ToOptimizeFunctions = true;
SRD.UseParallelizedSimplification = false;

% SymbolicEngine = SRD.GetSymbolicEngine;
% SymbolicEngine.LinearizationType = 'Naive';
% SRD.SaveSymbolicEngine(SymbolicEngine);

timerVal = tic;
SRD.DeriveEquationsForSimulation(ToLinearize, ToSimplify);
toc(timerVal);