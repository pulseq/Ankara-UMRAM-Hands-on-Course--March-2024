%% Pulseq tutorial for Ankara, 25.03.2024. Qingping Chen

%% Tasks:
% FID produces T2* weighted signal which is sensitive to macroscopic
% inhomogeneity (e.g. B0 inhomogeneity). ex02_fid2se is FID sequence adapted 
% from ex01_fid. Could you insert a 180-deg RF pulse (rf_ref) between 
% the excitation pulse (ref_ex) and ADC to generate a spin echo (SE) to reduce 
% the macroscopic inhomogeneity effect?

clear ; close all; clc ;
userID = 'chen' ; % your last name here

%% Define system properties
system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6) ;

%% Create a new sequence object
seq = mr.Sequence(system) ;

%% Set sequence parameters
Nrep = 1 ;% number of repetitions
Nx = 6400 ; % number of samples
adcDur = 128e-3 ; % ADC duration
rfDur = 1000e-6 ; % RF duration
TR = 500e-3 ; % unit: s
TE = 200e-3 ; % unit: s

%% Create a non-selective and refocusing pulse 
rf_ex = mr.makeBlockPulse(pi/2, 'Duration', rfDur, 'system', system) ;
rf_ref = mr.makeBlockPulse(pi,'Duration',rfDur, 'system', system, 'use', 'refocusing') ;

%% Define delays and ADC events
% recalculate delayTE1 and delayTE2 for rf_ref pulse
delayTE1 = TE/2 - rf_ex.shape_dur/2 - rf_ex.ringdownTime ;
delayTE2 = TE/2 ;
% delayTE1 = TE/2 - rf_ex.shape_dur/2 - rf_ex.ringdownTime - rf_ref.delay...
%     - rf_ref.shape_dur/2 ;
% delayTE2 = TE/2 - rf_ref.shape_dur/2 - rf_ref.ringdownTime - adcDur / 2 ;
adc = mr.makeAdc(Nx,'Duration',adcDur, 'system', system, 'delay', delayTE2) ;
delayTR = TR - mr.calcDuration(rf_ex) - delayTE1 - mr.calcDuration(adc) ;
% delayTR = TR - mr.calcDuration(rf_ex) - delayTE1 - mr.calcDuration(rf_ref) - mr.calcDuration(adc) ;
assert(delayTE1 >= 0) ;
assert(delayTE2 >= 0) ;
assert(delayTR >= 0) ;

%% Loop over repetitions and define sequence blocks
for i = 1:Nrep
    seq.addBlock(rf_ex) ;
    seq.addBlock(mr.makeDelay(delayTE1) ) ;
    % add rf_ref pulse to event block here
    seq.addBlock(adc) ;
    seq.addBlock(mr.makeDelay(delayTR) ) ;
end

%% Check whether the timing of the sequence is compatible with the scanner
[ok, error_report]=seq.checkTiming ;

if (ok)
    fprintf('Timing check passed successfully\n') ;
else
    fprintf('Timing check failed! Error listing follows:\n') ;
    fprintf([error_report{:}]) ;
    fprintf('\n') ;
end

seq.write(['ex02_fid2se_' userID '.seq']) ;       % Write to pulseq file

%% Plot sequence diagram
seq.plot() ;

%% plot k-spaces
% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP() ;
figure; plot(ktraj(1,:),ktraj(2,:),'b') ; % a 2D plot
axis('equal') ; % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.') ; % plot the sampling points
