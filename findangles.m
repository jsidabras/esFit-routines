function y = findangles(Sys,Exp,Opt);
    cori0 = [Sys.coria Sys.corib Sys.coric]*pi/180;          % initial crystal orientation
    nL = [1;0;0];              % rotation axis along lab x axis
    
    Sys.lwpp = 0.7;
    
    Sys1 = Sys;
    Sys1.gFrame = [Sys.gframea1 Sys.gframeb1 Sys.gframec1]*pi/180;
    Sys2 = Sys;
    Sys2.gFrame = [Sys.gframea1 Sys.gframeb1 Sys.gframec1]*pi/180 + [Sys.gframea2 Sys.gframeb2 Sys.gframec2]*pi/180;

    Exp.nPoints = 2048;
    % needed for total number of for loops
    files = dir('corrected_*.DTA');
    scalething = 2;
    lenme = length(files)/scalething;
    yq = [];
    cori1 = zeros(3,lenme);
    for i=1:lenme
        rho1 = (i-1)*5*scalething*pi/180 + Sys.pmang(i)*pi/180;
        cori1(:,i) = rotatecrystal(cori0,nL,rho1);
    end
    freq = [9.480429 9.480459 9.480461 9.480404 9.480403 9.480412 9.480417...
        9.480404 9.480395 9.480410 9.480410 9.480423 9.480423 9.480373...
        9.480384 9.480400 9.480402 9.479744 9.480625 9.480697 9.481300...
        9.481302 9.481298 9.481241 9.481993 9.481153 9.481152 9.481158...
        9.481160 9.481145 9.481157 9.481126 9.481129 9.481138 9.481142...
        9.481146];
    
    freq = freq(1:1*scalething:end);   
    parfor i=1:lenme
         PExp1 = struct('mwFreq', freq(i),'Range', Exp.Range,'CrystalSymmetry',Exp.CrystalSymmetry,'Harmonic', 0,'MolFrame', Sys.MolFrameA,'nPoints', Exp.nPoints,'CrystalOrientation', cori1(:,i));
         PExp2 = struct('mwFreq', freq(i),'Range', Exp.Range,'CrystalSymmetry',Exp.CrystalSymmetry,'Harmonic', 0,'MolFrame', Sys.MolFrameB,'nPoints', Exp.nPoints,'CrystalOrientation', cori1(:,i));

        [x,y11] = pepper(Sys1,PExp1,Opt);  
        [x,y12] = pepper(Sys2,PExp2,Opt);
                
        y1r(:,i) = rescale(y11 + y12, 'maxabs');    
    end
    for i=1:lenme
        yq = [yq(:); y1r(:,i)];
    end
    y = yq;
    
end