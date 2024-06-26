%% Pulseq tutorial for Ankara, 26.03.2024. Qingping Chen

%% Tasks:
% This exercise starts with ex15_gre2d_RFspoil. Could you try to adapt this
% sequence from 2D GRE to 2D radial?
% Hint: remove gyPre and gyReph. rotate gxPre, gx, and gxSpoil around
% z-axis:
% seq.addBlock(mr.rotate('z', phi, gxPre, gzReph)) ;
% seq.addBlock(mr.rotate('z', phi, gx, adc) ) ;
% seq.addBlock(mr.rotate('z', phi, gxSpoil, gzSpoil, mr.makeDelay(delayTR))) ;
clear ; close all; clc ;
userID = 'chen' ; % your last name here

%% set system limits and parameters
sys = mr.opts('MaxGrad', 22, 'GradUnit', 'mT/m', ...
    'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 20e-6) ;

seq = mr.Sequence(sys) ;           % Create a new sequence object
fov = 256e-3 ; Nx = 256 ; Ny = 256 ;     % Define FOV and resolution
alpha = 10 ;                       % flip angle
Nslices = 1 ;
sliceThickness = 3e-3 ;            % slice
sliceGap = 0 ;
TR = 20e-3 ;                       % repetition time TR
TE = 6e-3 ;                        % echo time TE
roDuration = 5.12e-3 ;             % ADC duration
delta = pi/Ny ;                    % angular increment; try golden angle pi*(3-5^0.5) or 0.5 of it

%% Create alpha-degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',3e-3,...
    'SliceThickness',sliceThickness,'apodization',0.42,'timeBwProduct',4,'system',sys) ;

%% Define other gradients and ADC events
deltak = 1/fov ;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',roDuration,'system',sys) ;
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys) ;
gxPre = mr.makeTrapezoid('x','Area',-gx.amplitude*(adc.dwell*(adc.numSamples/2+0.5)+0.5*gx.riseTime),...
    'Duration',1e-3,'system',sys) ;
gzReph = mr.makeTrapezoid('z','Area',-gz.area/2,'Duration',1e-3,'system',sys) ;
phaseAreas = ((0:Ny-1) - Ny/2) * deltak ;
gy = mr.makeTrapezoid('y','Area',max(abs(phaseAreas)),'Duration',mr.calcDuration(gxPre),'system',sys) ;
peScales = phaseAreas/gy.area ;
clear gyPre ;
for iY = 1:Ny
    gyPre(iY) = mr.scaleGrad(gy, peScales(iY)) ;
    gyReph(iY) = mr.scaleGrad(gy, -peScales(iY)) ;
end

% Prepare label for phase encoding
PElbl = 0:(Ny-1) ; % Caution: C-style numbering. starting from 0
% slice positions
slicePositions = (sliceThickness+sliceGap) * ((0:(Nslices-1)) - (Nslices-1)/2) ;

%% gradient spoiling
gxSpoil = mr.makeTrapezoid('x', 'Area', 2*Nx*deltak,'system',sys) ;
gzSpoil = mr.makeTrapezoid('z', 'Area', 4/sliceThickness,'system',sys) ;

%% Calculate timing
delayTE = ceil((TE - mr.calcDuration(gxPre) - gz.fallTime - gz.flatTime/2 ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime) * seq.gradRasterTime ;
delayTR = ceil((TR - mr.calcDuration(gz) - mr.calcDuration(gxPre) ...
    - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime ;
assert(delayTE >= 0 ) ;
assert(delayTR >= mr.calcDuration(gxSpoil, gyReph(1), gzSpoil) ) ;
% assert(delayTR >= mr.calcDuration(gxSpoil, gzSpoil) ) ;


%% Loop over phase encodes and define sequence blocks
tic ;
for i=1:Ny
    % RF spoiling (vary RF phase pseudo-randomly)
    rand_phase = mod(117*(i^2 + i + 2), 360) * pi/180 ;
    rf.phaseOffset = rand_phase ;
    adc.phaseOffset = rand_phase ;
    seq.addBlock(rf, gz) ; % slice-selective excitation
    phi = delta * (i-1) ;
    seq.addBlock(gxPre, gyPre(i), gzReph, mr.makeLabel('SET', 'LIN', PElbl(i))) ;
    seq.addBlock(mr.makeDelay(delayTE)) ;
    seq.addBlock(gx, adc) ;
    seq.addBlock(mr.makeDelay(delayTR), gxSpoil, gyReph(i), gzSpoil) ;
end
toc ;

%% check whether the timing of the sequence is correct
[ok, error_report] = seq.checkTiming ;

if (ok)
    fprintf('Timing check passed successfully\n') ;
else
    fprintf('Timing check failed! Error listing follows:\n') ;
    fprintf([error_report{:}]) ;
    fprintf('\n') ;
end

%% prepare sequence export
seq.setDefinition('Name', 'gre') ;
% the following definitions have effect in conjunction with LABELs 
seq.setDefinition('FOV', [fov fov max(slicePositions)-min(slicePositions)+sliceThickness]) ;
seq.setDefinition('SlicePositions', slicePositions) ;
seq.setDefinition('SliceThickness', sliceThickness) ;
seq.setDefinition('SliceGap', sliceGap) ;

seq.write(['ex24_gre2radial_' userID '.seq']) ;       % Write to pulseq file

%% evaluate label settings
% adc_lbl = seq.evalLabels('evolution','adc') ;
% figure ;
% hold on ;
% plot(adc_lbl.LIN) ;
% legend('lin') ;
% title('evolution of labels/counters') ;

%% plot sequence and k-space diagrams

seq.plot('timeRange', [0 5]*TR) ;

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP() ;

% plot k-spaces
figure; plot(ktraj(1,:),ktraj(2,:),'b') ; % a 2D plot
axis('equal') ; % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.') ; % plot the sampling points



