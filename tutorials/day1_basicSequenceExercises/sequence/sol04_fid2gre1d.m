%% Pulseq tutorial for Ankara, 25.03.2024. Qingping Chen

%% Tasks:
% Gradient enables spatial encoding and thus multi-dimensional images can 
% be produced. We can also use gradient for dephasing and recall an echo, 
% aka gradient-recalled echo (GRE). This exercise starts with FID with
% slice-selective excitation. Could you try to insert a readout gradient
% and its prephaser to generate a echo?

clear ; close all; clc ;
userID = 'chen' ; % your last name here

%% set system limits and parameters
sys = mr.opts('MaxGrad', 22, 'GradUnit', 'mT/m', ...
    'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 20e-6) ;

seq = mr.Sequence(sys) ;           % Create a new sequence object
fov = 256e-3 ; Nx = 256 ; Ny = 1 ;     % Define FOV and resolution
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

%% Calculate timing
delayTE = ceil((TE - mr.calcDuration(gxPre) - gz.fallTime - gz.flatTime/2 ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime) * seq.gradRasterTime ;
delayTR = ceil((TR - mr.calcDuration(gz) - mr.calcDuration(gxPre) ...
    - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime ;
assert(delayTE >= 0 ) ;
assert(delayTR >= 0 ) ;

%% Loop over phase encodes and define sequence blocks
tic ;
for i=1:Ny
    seq.addBlock(rf, gz) ; % slice-selective excitation
    seq.addBlock(gxPre, gzReph) ;
    seq.addBlock(mr.makeDelay(delayTE)) ;
    seq.addBlock(gx, adc) ;
    seq.addBlock(mr.makeDelay(delayTR)) ;
end
toc ;

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming ;

if (ok)
    fprintf('Timing check passed successfully\n') ;
else
    fprintf('Timing check failed! Error listing follows:\n') ;
    fprintf([error_report{:}]) ;
    fprintf('\n') ;
end

%% prepare sequence export
seq.setDefinition('FOV', [fov fov sliceThickness]) ;
seq.setDefinition('Name', 'gre1d') ;

seq.write(['sol04_fid2gre1d_' userID '.seq']) ;       % Write to pulseq file

%% plot sequence and k-space diagrams

seq.plot('timeRange', [0 5]*TR) ;

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP() ;

% plot k-spaces
figure; plot(ktraj(1,:),ktraj(2,:),'b') ; % a 2D plot
axis('equal') ; % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.') ; % plot the sampling points



