% spectra from single crystal rotation
%===========================================

clear, clf

% Spin parameters

% Experimental parameters
Exp.mwFreq = 9.48075;
Exp.Range = [312.5 362.4756];
Exp.CrystalSymmetry = 'P1211';
Exp.Harmonic = 0;

Sys.g = [2.1 2.04 1.99];

vals = [-9.9905 -104.3957 179.1777 -143.6922 -89.5230 139.4536 -5.1538 -6.9590 -3.3403]; % 0.0652
 BestVals = struct('coria', vals(1),'corib', vals(2),'coric', vals(3),...
      'gframea1', vals(4),'gframeb1', vals(5),'gframec1', vals(6),...
      'gframea2', vals(7),'gframeb2', vals(8),'gframec2', vals(9));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global Fitting
% Use Corrected data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spc_cont = [];
% files = dir('corrected_*.DTA');
% for i=1:length(files)
%     [B,spc1,Params] = eprload(files(i).name);
%     spc1r = rescale(real(spc1), 'maxabs');
%     spc_cont = [spc_cont(:); spc1r(:)];
% end
% 
% spc_cont = rescale(spc_cont, 'maxabs');
% Use gaussians of known width at measured peaks
spc_cont = csvread('PerfectData.csv');

% MolFrame calculated from PCB ID: 4XDC
R1=[-0.3308 -0.7701 +0.5454;
    -0.9376 +0.2027 -0.2825;
    +0.1071 -0.6049 -0.7891];
    
R2=[-0.5948 +0.4557 -0.6622;
    -0.6664 -0.7403 +0.0891;
    -0.4496 +0.4943 +0.7440];
Sys.MolFrameA=eulang(R1);
Sys.MolFrameB=eulang(R2);

% -CN from Adamska-Venkatesh 10.1039/C4CP05426A
AN6 = [-1.3 -1.1 6.2]; N6_AFrame = [0 50 90]*pi/180;
Sys = nucspinadd(Sys,'14N',AN6,N6_AFrame);

Sys.coria = BestVals.coria;
Sys.corib = BestVals.corib;
Sys.coric = BestVals.coric;
Vary.coria = 90;
Vary.corib = 90;
Vary.coric = 90;

Sys.gframea1 = BestVals.gframea1;
Sys.gframeb1 = BestVals.gframeb1;
Sys.gframec1 = BestVals.gframec1;
Vary.gframea1 = 90;
Vary.gframeb1 = 90;
Vary.gframec1 = 90;

Sys.gframea2 =  BestVals.gframea2;
Sys.gframeb2 =  BestVals.gframeb2;
Sys.gframec2 =  BestVals.gframec2;
Vary.gframea2 = 7.5;
Vary.gframeb2 = 7.5;
Vary.gframec2 = 7.5;

Sys.pmang = zeros(1,18)
% Vary.pmang = ones(1,18)*3.5

FitOpt.Scaling = 'lsq';
FitOpt.Method = 'swarm';     % for the genetic algorithm
FitOpt.PopulationSize = 1000;
FitOpt.maxGenerations = 500;
FitOpt.nParticles = 50000;
FitOpt.RandomStart = 1;
FitOpt.nTrials = 100;
% FitOpt.maxTime = 90;
% [BestVals, BestFit] = esfit('findangles',spc_cont,Sys,Vary,Exp,[],FitOpt);
esfit('findangles',spc_cont,Sys,Vary,Exp,[],FitOpt)
