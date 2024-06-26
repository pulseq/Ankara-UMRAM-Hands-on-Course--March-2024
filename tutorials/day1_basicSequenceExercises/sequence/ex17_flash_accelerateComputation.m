%% Pulseq tutorial for Ankara, 25.03.2024. Qingping Chen

%% Tasks:
% This exercise starts with ex16_flash_fasterTiming. The computation time
% for building a sequence might be very long in case that the sequence is
% complicated. this is because the seq.addBlock() is time-consuming for
% searching the phase/magnitude shape in the library. Try to accelerate teh
% computation time by register the unchanged gradient/RF event to the
% library in advance. Hint: seq.registerGradEvent(), seq.registerRfEvent().
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
roDuration = 5.12e-3 ;              % ADC duration

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
% gxSpoil = mr.makeTrapezoid('x', 'Area', 2*Nx*deltak,'system',sys) ;
% we cut the RO gradient into two parts for the optimal spoiler timing
[gx1, ~] = mr.splitGradientAt(gx, gx.riseTime + gx.flatTime) ;

% gradient spoiling
gxSpoil = mr.makeExtendedTrapezoidArea(gx.channel, gx.amplitude, 0, 2*Nx*deltak, sys) ;
gxSpoil.delay = mr.calcDuration(gx1) ;
gx_add = mr.addGradients({gx1, gxSpoil},'system',sys) ;

gzSpoil = mr.makeTrapezoid('z', 'Area', 4/sliceThickness,'system',sys) ;

% add delay of gx1 duration to gzSpoil and gyReph(i)
gzSpoil.delay = max(mr.calcDuration(gx1), gzSpoil.delay) ;
for iY = 1:Ny
    gyReph(iY).delay = max(mr.calcDuration(gx1), gyReph(iY).delay) ;
end

%% Calculate timing
delayTE = ceil((TE - mr.calcDuration(gxPre) - gz.fallTime - gz.flatTime/2 ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime) * seq.gradRasterTime ;
% delayTR = ceil((TR - mr.calcDuration(gz) - mr.calcDuration(gxPre) ...
%     - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime ;
assert(delayTE >= 0 ) ;
% assert(delayTR >= mr.calcDuration(gxSpoil, gyReph(1), gzSpoil) ) ;

%% accelerate computations
% preregister constant objects to accelerate computations
% this is not necessary, but accelerates the sequence creation
gxPre.id = seq.registerGradEvent(gxPre) ;
gx_add.id = seq.registerGradEvent(gx_add) ;
gzSpoil.id = seq.registerGradEvent(gzSpoil) ;
% the phase of the RF object will change, therefore we only per-register the shapes
rf.id = seq.registerRfEvent(rf) ; % task: try to correct this line?

for iY = 1:Ny
    gyPre(iY).id = seq.registerGradEvent(gyPre(iY)) ;
    gyReph(iY).id = seq.registerGradEvent(gyReph(iY)) ;
end

%% Loop over phase encodes and define sequence blocks
tic ;
for i=1:Ny
    % RF spoiling (vary RF phase pseudo-randomly)
    rand_phase = mod(117*(i^2 + i + 2), 360) * pi/180 ;
    rf.phaseOffset = rand_phase ;
    adc.phaseOffset = rand_phase ;
    seq.addBlock(rf, gz) ; % slice-selective excitation
    seq.addBlock(gxPre, gyPre(i), gzReph, mr.makeLabel('SET', 'LIN', PElbl(i))) ;
    seq.addBlock(mr.makeDelay(delayTE)) ;
    seq.addBlock(gx_add, adc, gyReph(i), gzSpoil) ;
%     seq.addBlock(mr.makeDelay(delayTR), gxSpoil, gyReph(i), gzSpoil) ;
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

seq.write(['ex17_flash_accelerateComputation_' userID '.seq']) ;       % Write to pulseq file

%% evaluate label settings
adc_lbl = seq.evalLabels('evolution','adc') ;
figure ;
hold on ;
plot(adc_lbl.LIN) ;
legend('lin') ;
title('evolution of labels/counters') ;

%% plot sequence and k-space diagrams

seq.plot('timeRange', [0 5]*TR) ;

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP() ;

% plot k-spaces
figure; plot(ktraj(1,:),ktraj(2,:),'b') ; % a 2D plot
axis('equal') ; % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.') ; % plot the sampling points