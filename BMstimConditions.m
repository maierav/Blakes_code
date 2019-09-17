%% mainframe rewrite
%%script to select an input file amongst a decision tree, load in .nev,
%%.ns2, and ns6, pull out stim onset and tie those timepoints to the
%%raw neural data. Each stim onset with animal fixation is a trial. Generate LFP, aMUA, and CSD --triggered to a reference
%%window. Average across trials. Baseline correct the data. Plot each
%%contact channel from -50 to 300ms relative to stim onset. 
clear
%% Establish directories and set path

if strcmp(getenv('USER'),'maierav')                                      %retrieves environment variable 'USER' 
    npmkdir  = '/Users/alex 1/Desktop/LAB/Brock/OLD/NPMK-4.5.3.0/NPMK/'; %directory for Alex's machine
    nbanalysisdir   = '/Users/alex 1/Desktop/LAB/bootcamp/nbanalysis/';  %directory for Alex's machine
    datadir  = '/Users/alex 1/Desktop/LAB/';                             %directory for the stored data
else
    npmkdir  = '/users/bmitc/Documents/MATLAB/NPMK/';                    %neural processing matlab kit (NPMK)
    nbanalysisdir   = '/users/bmitc/Documents/MATLAB/nbanalysis/';       %directory with various tools for opening, loading, and processing 
    datadir  = '/users/bmitc/Documents/MATLAB/data/';
end

addpath(genpath(npmkdir))
addpath(genpath(nbanalysisdir))
addpath(genpath(datadir))

BRdatafile = '161006_E_rfori001'; 
filename = [datadir BRdatafile];

%% Define stimulus patterns and select from among them

patterns   = {'rforidrft','rfsfdrft','posdisparitydrft','disparitydrft','cinterocdrft','coneinterocdrft','conedrft', ...
                'colorflicker','bwflicker','rfori','rfsize','cinteroc','color','rfsf','mcosinteroc','dotmapping'}; 
for p = 1:length(patterns) 
    
    pattern = patterns{p};    %pattern is individual elements of the above array
    
    if any(strfind(BRdatafile,pattern))         %if BRdatafile contains any string the same as pattern
        startlog = strfind(BRdatafile,pattern); %create a variable to store the number of those patches
        if ~isequal(BRdatafile(startlog:end-3),pattern), continue
        else
        match = patterns{p};
        end
    end
    
end

if isequal(match,'dotmapping')
    ext = '.gDotsXY_di';
else
    ext = ['.g' upper(match) 'Grating_di'];
end

if contains(ext,'RFORI')
    grating = readgGrating([filename ext]);
else
    error('add more extensions to choose from')
end

%% Load event times and codes

NEV             = openNEV([filename '.nev'],'noread','overwrite');
EventCodes      = NEV.Data.SerialDigitalIO.UnparsedData - 128;          %we don't know why we subtract 128
EventSamples    = NEV.Data.SerialDigitalIO.TimeStamp;                   %Events in samples 
EventTimes      = floor(NEV.Data.SerialDigitalIO.TimeStampSec.*1000);   %floor rounds to nearest integer and then convert event to ms 
[pEvC, pEvT]    = parsEventCodesML(EventCodes,EventSamples);            %sorts codes, samps or times into trials

%The following is a structure that only contains trials where the animal
%did not break fixation, and includes grating info and stimulus onsets. 

STIM            = sortStimandTimeData(grating,pEvC,pEvT,'stim'); 
STIM.onsetsdown = floor(STIM.onsets./30);

%% Load LFP with NS2 file
clear ext
ext = 'ns2';

%first retrieve sort direction and bank info
NS_Header    = openNSx(strcat(filename,'.',ext),'noread');
banks        = unique({NS_Header.ElectrodesInfo.ConnectorBank}); banks(ismember(banks,'E')) = []; % bank E is BNC cable inputs

for b = 1:length(banks)
    clear neural label 
    neural       = strcmp({NS_Header.ElectrodesInfo.ConnectorBank},banks{b}); 
    firstlabel   = cell2mat({NS_Header.ElectrodesInfo(find(neural,1,'first')).Label}); 
    if str2double(firstlabel(3:4)) < 2
        sortdirection = 'ascending'; 
    else
        sortdirection = 'descending'; 
    end
end

if any(strfind(firstlabel,'eD'))
    ebank = 'eD'; %electrode bank
else 
    ebank = 'eB';
end

% load LFP with NS2 file

lfp = getLFP(filename,ext,ebank,sortdirection);

%% Create a new STIM.field for separate orientation tilts - work in progress
clear i count1 count2
count1 = 0; count2 = 0;
for i = 1:size(STIM.tilt,1)
    if STIM.tilt(i) == 157.5 && STIM.eye(i) == 3
        count1 = count1+1;
        STIM.ori1onsets(count1,:) = STIM.onsetsdown(i);
    elseif STIM.tilt(i) == 0 && STIM.eye(i) == 3
        count2 = count2+1;
        STIM.ori2onsets(count2,:) = STIM.onsetsdown(i);
    end    
end

%% LOAD LFP and analog MUA with NS6 file.

clear ext NS_header banks neural 
% Read in NS Header
ext          = 'ns6'; 
NS_Header    = openNSx(strcat(filename,'.',ext),'noread');

% get basic info about recorded data
neural       = strcmp({NS_Header.ElectrodesInfo.ConnectorBank},ebank(2)); % logicals where contact bank name matches electrode of interest
N.neural     = sum(neural); % number of neural channels 
NeuralLabels = {NS_Header.ElectrodesInfo(neural).Label}; %get labels
Fs           = NS_Header.MetaTags.SamplingFreq; % get sampling frequency
nyq          = Fs/2; 
r            = Fs/1000; 

% counters
clear nct
nct = 0;

% process data electrode by electrode
for e = 1:length(neural)
    
    if neural(e) == 1    % neural is a vector of logicals, 1 = contacts we want

        nct = nct+1;
        
        % open data for this channel. 
        clear NS DAT
        electrode = sprintf('c:%u',e);
        NS        = openNSx(strcat(filename,'.',ext),electrode,'read','uV');
        DAT       = NS.Data; NS.Data = [];  % this is the whole signal on one channel, 30 kHz!
        
        
        % preallocate data matrices 
        if nct == 1
            N.samples = length(DAT); 
            LFP       = nan(ceil(N.samples/r),N.neural); % preallocating for downsampled data
            MUA       = nan(ceil(N.samples/r),N.neural);
        end
        
        % extract the LFP. 
        clear lpc lWn bwb bwa lpLFP
        lpc       = 200; %low pass cutoff
        lWn       = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low');
        lpLFP      = filtfilt(bwb,bwa,DAT);  %low pass filter 
        
        % extract the MUA:
        clear hpc hWn bwb bwa hpMUA
        hpc       = 750;  %high pass cutoff
        hWn       = hpc/nyq;
        [bwb,bwa] = butter(4,hWn,'high');
        hpMUA     = filtfilt(bwb,bwa,DAT); %high pass filter
        
        % low pass at 5000 Hz and rectify (take the absolute value) 
        clear lpc lWn bwb bwa 
        lpc       = 5000;  % cutoff
        lWn       = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low');
        hpMUA     = abs(filtfilt(bwb,bwa,hpMUA)); %low pass filter &rectify
        
        % low pass filter at x Hz. 
        clear lpc lWn bwb bwa lpMUA
        lpc       = 200; %low pass cutoff
        lWn       = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low'); 
        lpMUA     = filtfilt(bwb,bwa,hpMUA);  %low pass filter to smooth
        
      
        % decimate both LFP and analog MUA (aMUA) to get 1kHz samp freq
        MUA(:,nct) = decimate(lpMUA,r); 
        LFP(:,nct) = decimate(lpLFP,r); 
        
        clear DAT 
        
    end
    
end

% sort electrode contacts in ascending order:
idx    = nan(length(NeuralLabels),1); 
for ch = 1:length(NeuralLabels)
    chname  = strcat(sprintf('%s',ebank),sprintf('%02d',ch));
    idx(ch) = find(contains(NeuralLabels,chname));
    
end

switch sortdirection
    case 'ascending'
        MUA = MUA(:,idx);
        LFP = LFP(:,idx);
        sortedLabels = NeuralLabels(idx); 
    case 'descending'
        MUA = MUA(:,fliplr(idx));
        LFP = LFP(:,fliplr(idx));
        sortedLabels = NeuralLabels(fliplr(idx)); 
end

%% calculate CSD 
% calculate CSD before triggering to trials OR on the trial data BUT not on
% the mean LFP. 

CSD = mod_iCSD(LFP')';  % this function takes LFP in channels x samples so let's transpose LFP and then flip it right back 
                        % feed in units of microV and get back units of
                        % nA/mm^3


CSD = padarray(CSD,[0 1],NaN,'replicate'); % pad array if you want to keep the matrix the same size on the channel
                                           % dimension as the other matrices

%% trigger the neural data to the event codes of interest                                         

pre   = -50;
post  = 300; 

STIM.LFP_ori1  = trigData(LFP,STIM.ori1onsets,-pre,post); %pre variable is in absolute units 
STIM.CSD_ori1  = trigData(CSD,STIM.ori1onsets,-pre,post); 
STIM.aMUA_ori1 = trigData(MUA,STIM.ori1onsets,-pre,post);

STIM.LFP_ori2  = trigData(LFP,STIM.ori2onsets,-pre,post); 
STIM.CSD_ori2  = trigData(CSD,STIM.ori2onsets,-pre,post); 
STIM.aMUA_ori2 = trigData(MUA,STIM.ori2onsets,-pre,post); 

%% Averaging across trials
fields.STIM = fieldnames(STIM);
for avtr=1:numel(fields.STIM)
    AVG.(fields.STIM{avtr})  = mean(STIM.(fields.STIM{avtr}),3);
end

%% Baseline correction ERROR
% bl AVG
%fields.AVG = fieldnames(AVG);
%for blavg=1:numel(fields.AVG:end-6)
   % bl  = mean(AVG.(fields.AVG{blavg})(:,25:75),2);
   % BLavg.(fields.AVG{blavg})  = AVG.(fields.AVG{blavg}) - bl;
%end

%% Interpolate CSD

%% Plot! 
% LFP comparison; Ori 1 vs Ori 2

[bAVG_LFP_ori1] = BMbasecorrect(AVG.LFP_ori1);

[bAVG_LFP_ori2] = BMBMbasecorrect(AVG.LFP_ori2);

refwin = pre:post;

figure;
subplot(1,2,1)
f_ShadedLinePlotbyDepth(bAVG_LFP_ori1,1:24,refwin,[],1)
hold on
plot([0 0], ylim,'k')
plot([800 800], ylim,'k')
title({'LFP - 157.5 tilt',BRdatafile}, 'Interpreter', 'none')
xlabel('time (ms)')
ylabel('Contacts indexed down from surface')
hold off

subplot(1,2,2)
f_ShadedLinePlotbyDepth(bAVG_LFP_ori2,1:24,refwin,[],1)
hold on
plot([0 0], ylim,'k')
plot([800 800], ylim,'k')
title({'LFP - 0 tilt',BRdatafile}, 'Interpreter', 'none')
xlabel('time (ms)')
ylabel('Contacts indexed down from surface')
hold off

%% MUA comparison

[bAVG_aMUA_ori1] = BMbasecorrect(AVG.aMUA_ori1);

[bAVG_aMUA_ori2] = BMbasecorrect(AVG.aMUA_ori2);

figure;
subplot(1,2,1)
f_ShadedLinePlotbyDepth(bAVG_aMUA_ori1,1:24,refwin,[],1)
hold on
plot([0 0], ylim,'k')
plot([800 800], ylim,'k')
title({'aMUA - 157.5 tilt',BRdatafile}, 'Interpreter', 'none')
xlabel('time (ms)')
ylabel('Contacts indexed down from surface')
hold off

subplot(1,2,2)
f_ShadedLinePlotbyDepth(bAVG_aMUA_ori2,1:24,refwin,[],1)
hold on
plot([0 0], ylim,'k')
plot([800 800], ylim,'k')
title({'aMUA - 0 tilt',BRdatafile}, 'Interpreter', 'none')
xlabel('time (ms)')
ylabel('Contacts indexed down from surface')
hold off

%% CSD comparison

[bAVG_CSD_ori1] = BMbasecorrect(AVG.CSD_ori1);

[bAVG_CSD_ori2] = BMbasecorrect(AVG.CSD_ori2);

f_bAVG_CSD_ori1 = filterCSD(bAVG_CSD_ori1')';
f_bAVG_CSD_ori2 = filterCSD(bAVG_CSD_ori2')';

figure;
subplot(1,2,1)
imagesc(refwin,1:24,f_bAVG_CSD_ori1');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; v = vline(0); set(v,'color','k','linestyle','-.','linewidth',1);
set(gca,'tickdir','out');  
climit = max(abs(get(gca,'CLim'))*1);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
%plot([0 0], ylim,'k')
%plot([-800 800], ylim,'k')
title({'Interpolated CSD - 157.5 tilt',BRdatafile}, 'Interpreter', 'none')
xlabel('\fontsize{12}Time (ms)')
clrbar = colorbar; clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',12,'VerticalAlignment','middle');
ylabel('\fontsize{12}Contacts indexed down from the surface');
hold off

subplot(1,2,2)
imagesc(refwin,1:24,f_bAVG_CSD_ori2');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; v = vline(0); set(v,'color','k','linestyle','-.','linewidth',1);
set(gca,'tickdir','out');  
climit = max(abs(get(gca,'CLim'))*1);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
hold on;
%plot([0 0], ylim,'k')
%plot([-800 800], ylim,'k')
title({'Interpolated CSD - 0 tilt',BRdatafile}, 'Interpreter', 'none')
xlabel('\fontsize{12}Time (ms)')
clrbar = colorbar; clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',12,'VerticalAlignment','middle');
ylabel('\fontsize{12}Contacts indexed down from the surface');