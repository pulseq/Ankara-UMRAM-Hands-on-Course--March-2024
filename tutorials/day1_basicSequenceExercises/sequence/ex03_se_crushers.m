%% Pulseq tutorial for Ankara, 25.03.2024. Qingping Chen

%% Tasks:
% The 180-deg RF pulse is typically not perfect. A pair of crusher
% gradients can be used to suppress unwanted FID. In ex03_se_crushers,
% could you try to modify the spin-echo FID (ex02_se) to insert crushers
% before and after the 180-deg RF pulse?

clear ; close all; clc ;
userID = 'chen' ; % your last name here

%% Define system properties
system = mr.opts('MaxGrad',22,'GradUnit','mT/m',...
    		'MaxSlew',120,'SlewUnit','T/m/s',...
    		'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
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
spA = 1000; % spoiler area in 1/m (=Hz/m*s)

%% Create a non-selective and refocusing pulse 
rf_ex = mr.makeBlockPulse(pi/2, 'Duration', rfDur, 'system', system) ;
rf_ref = mr.makeBlockPulse(pi,'Duration',rfDur, 'system', system, 'use', 'refocusing') ;

%% calculate spoiler gradient
g_sp = mr.makeTrapezoid('z','Area',spA,'system',system) ;
% rf_ref.delay = max(mr.calcDuration(g_sp), rf_ref.delay) ;

%% Define delays and ADC events
delayTE1 = TE/2 - rf_ex.shape_dur/2 - rf_ex.ringdownTime - rf_ref.delay...
    - rf_ref.shape_dur/2 ;
delayTE2 = TE/2 - rf_ref.shape_dur/2 - rf_ref.ringdownTime - adcDur / 2 ;
adc = mr.makeAdc(Nx,'Duration',adcDur, 'system', system, 'delay', delayTE2) ;
delayTR = TR - mr.calcDuration(rf_ex) - delayTE1 - mr.calcDuration(rf_ref) - mr.calcDuration(adc) ;
assert(delayTE1 >= 0) ;
assert(delayTE2 >= mr.calcDuration(g_sp)) ; %% ADC delay > g_sp
assert(delayTR >= 0) ;

%% Loop over repetitions and define sequence blocks
for i = 1:Nrep
    seq.addBlock(rf_ex) ;
    seq.addBlock(mr.makeDelay(delayTE1) ) ;
    seq.addBlock(rf_ref) ;
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

seq.write(['ex03_se_crushers_' userID '.seq']) ;       % Write to pulseq file


%% Plot sequence diagram
seq.plot() ;

%% plot k-spaces
% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP() ;
figure; plot(ktraj(1,:),ktraj(2,:),'b') ; % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.') ; % plot the sampling points

