% spectra from single crystal rotation
%===========================================
clear, clf
% vals = [-9.8897 -97.2763 177.0524 -121.2711 -106.2582 154.1617 -3.9470 -4.7960 0.0208]; % 0.0754 RMSD
% vals = [-9.4574 -104.6001 179.2189 -145.3516 -88.3269 142.4581 -2.6790 -9.0206 -8.6355]; % 0.0729 RMSD
vals = [-9.9905 -104.3957 179.1777 -143.6922 -89.5230 139.4536 -5.1538 -6.9590 -3.3403]; % 0.0652 RMSD
BestVals = struct('coria', vals(1),'corib', vals(2),'coric', vals(3), 'gframea1', vals(4),'gframeb1', vals(5),'gframec1', vals(6), 'gframea2', vals(7),'gframeb2', vals(8),'gframec2', vals(9));
  
% Experimental parameters
Exp.mwFreq = 9.48075;
Exp.Range = [312.5 362.4756];
Exp.CrystalSymmetry = 'P1211'; % from PDB ID: 4XDC
Exp.Harmonic = 0;


% Spin parameters
Sys.g = [2.1 2.04 1.99];
Sys.lwpp = [0.7];

% Both MolFrames are calculated from PDB ID: 4XDC
R1=[-0.3308 -0.7701 +0.5454;
    -0.9376 +0.2027 -0.2825;
    +0.1071 -0.6049 -0.7891];
R2=[-0.5948 +0.4557 -0.6622;
    -0.6664 -0.7403 +0.0891;
    -0.4496 +0.4943 +0.7440];
Sys.MolFrameA=eulang(R1);
Sys.MolFrameB=eulang(R2);

% Generate orientations in a single rotation plane
nL = [1; 0; 0;];
cori0 = [BestVals.coria BestVals.corib BestVals.coric]*pi/180;
rho1 = [0:5:180]*pi/180;
cori1 = rotatecrystal(cori0,nL,rho1);
Exp.CrystalOrientation = cori1;

% Duplicate Exp split for the 2 assymetric units.
Exp1=Exp;
Exp2=Exp;
Exp1.MolFrame = Sys.MolFrameA;
Exp2.MolFrame = Sys.MolFrameB;

% Duplicate Sys split for the 2 assymetric units.
Sys1 = Sys;
Sys1.gFrame = [BestVals.gframea1 BestVals.gframeb1 BestVals.gframec1]*pi/180;
Sys2 = Sys;
% Assume gFrameA is related to gFrameB by: gFrameA + delta-angle
Sys2.gFrame = [BestVals.gframea1 BestVals.gframeb1 BestVals.gframec1]*pi/180 + [BestVals.gframea2 BestVals.gframeb2 BestVals.gframec2]*pi/180;

% Simulate roadmape
Opt.Output = 'separate';  % make sure spectra are not added up
Bres1 = resfields(Sys1,Exp1,Opt);
Bres2 = resfields(Sys2,Exp2,Opt);

Y1 = [0:5:180]+2;
X1 = Bres1(1,:);
X2 = Bres1(2,:);
X3 = Bres2(1,:);
X4 = Bres2(2,:);

% plotting
figure1 = figure;

YMatrix1 = [];
files = dir('corrected_*.DTA');
for i=1:length(files)
    [B,spc1,Params] = eprload(files(i).name);
    spc1r = rescale(real(spc1), 'maxabs')*5 + 5*(i-1);
    YMatrix1 = [YMatrix1 spc1r];
    X5 = B/10;
end


% Create axes
axes1 = axes('Parent',figure1,'ColorOrder',[0 0 0]);
hold(axes1,'on');

% Create plot
plot(X1,Y1,'LineWidth',3,'LineStyle','--',...
    'Color',[1 0.08 0.65]);

% Create plot
plot(X2,Y1,'LineWidth',3,'LineStyle','--',...
    'Color',[0.47 0.67 0.19]);

% Create plot
plot(X3,Y1,'LineWidth',3,'LineStyle','--','Color',[0 0 1]);

% Create plot
plot(X4,Y1,'LineWidth',3,'LineStyle','--','Color',[1 0 0]);

% Create multiple lines using matrix input to plot
plot(X5,YMatrix1,'LineWidth',2,'Color',[0 0 0],'Parent',axes1);

% Create xlabel
xlabel({'Magnetic Field [mT]'});
ylabel({'Theta [deg]'});
hold off
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[328 342]);
ylim(axes1,[-5 185]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1, 'FontSize', 24, 'GridAlpha', 0.25, 'GridLineStyle', '--', 'LineStyleOrder', {'-'}, 'LineWidth', 2, 'TickDir', 'out', 'XColor', [0 0 0], 'XGrid', 'on', 'XMinorTick', 'on', 'YColor', [0 0 0],'DataAspectRatio',[1 10 1]);