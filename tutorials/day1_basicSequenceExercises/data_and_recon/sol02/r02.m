clear all; close all; clc;
% very basic FID data handling example
%
% needs mapVBVD in the path

%% Load data from the given directory sorted by name
path = pwd ; % directory to be scanned for data files
nF = 1 ; % the number of the file / data set to load

pattern = '*meas*.dat' ;
D = dir([path filesep pattern]) ;
[~,I] = sort(string({D(:).name})) ;
data_file_path = [path filesep D(I(nF)).name] ;
disp(['loading `' data_file_path '´ ...']) ;
twix_obj = mapVBVD(data_file_path) ;

%% raw data preparation
if iscell(twix_obj)
    data_unsorted = double(twix_obj{end}.image.unsorted()) ;
else
    data_unsorted = double(twix_obj.image.unsorted()) ;
end

channels = size(data_unsorted, 2) ;
adc_len = size(data_unsorted, 1) ;
readouts = size(data_unsorted, 3) ;

rawdata = permute(data_unsorted, [1,3,2]) ;

%% Load sequence from the file with the same name as the raw data
seq_file_path = 'sol02_fid2se_chen.seq' ;
fprintf(['loading `' seq_file_path '´ ...\n']) ;
seq = mr.Sequence() ;              % Create a new sequence object
if nF~=2 % 'detectRFuse' will not fork for multi-FA FID sequence because alpha excitations with the flip-angle larger than 90 degrees will be detected as refocusing pulses
    seq.read(seq_file_path,'detectRFuse') ;
else
    seq.read(seq_file_path) ;
end

%% we just want t_adc
[ktraj_adc, t_adc, ktraj, t, t_excitation, t_refocusing] = seq.calculateKspacePP() ;

t_adc = reshape(t_adc,[adc_len,readouts]) ;
% calc TE but account for dummy scans...
%[~,iND]=max(t_excitation(t_excitation<t_adc(1,1)));
%t_e=t_adc-t_excitation(ones(1,adc_len),iND:end);
t_relevant_excitation = zeros(1,readouts) ;
for i = 1:readouts
    t_relevant_excitation(i) = max(t_excitation(t_excitation<t_adc(1,i))) ;
end
t_e = t_adc-t_relevant_excitation(ones(1,adc_len),:) ;

%% plot raw data
figure; plot(t_e, abs(rawdata(:,:,1))); title('k-space data');
hold on;
if channels>1
    for j=2:channels
        plot(t_e, abs(rawdata(:,:,j))) ;
    end
end
 
plot(t_e(1), 0) ; % trick to force Y axis scaling
xlabel('time since excitation /s');
exportgraphics(gcf,'sol02_kspace.png',...
    'ContentType','vector',...
    'BackgroundColor','none') ;

if (readouts>1)
    figure; plot(abs(rawdata(4,:,1))); title('time evolution (4th FID point)');
    if channels>1
        hold on;
        for j=2:channels
            plot(abs(rawdata(4,:,j)));
        end
    end

    dataavg=squeeze(mean(rawdata,2));
    figure; plot(t_e, abs(dataavg(:,1))); title('averaged FIDs');
    if channels>1
        hold on;
        for j=2:channels
            plot(t_e, abs(dataavg(:,2)));
        end
    end
    xlabel('time since excitation /s');
end

%% plot spectr(um/a)
data_fft = fftshift(fft(fftshift(rawdata,1)), 1) ;
w_axis = (-adc_len/2:(adc_len/2-1))/((t_adc(end,1,1)-t_adc(1,1,1))/(adc_len-1)*adc_len)/twix_obj.hdr.Dicom.lFrequency*1e6' ;
figure ;
plot(w_axis, abs(data_fft(:,:,1)), 'lineWidth', 1) ;
title('abs spectr(um/a)') ;
xlim([-10 10]) ; xlabel('frequency /ppm') ;
set(gca,'Xdir','reverse') ;

exportgraphics(gcf,'sol02_im.png',...
    'ContentType','vector',...
    'BackgroundColor','none') ;