% working code for handling and analzying neural data
clear 

%% establish directories for toolkits and data
npmkdir    = '/users/bmitc/Documents/MATLAB/NPMK/'; 
nbanadir   = '/users/bmitc/Documents/MATLAB/nbanalysis/'; 

directory  = '/users/bmitc/Documents/MATLAB/data/';
BRdatafile = '161006_E_rfori001';       %sample data 
filename   = [directory BRdatafile]; 

addpath(genpath(directory))
addpath(genpath(npmkdir))
addpath(genpath(nbanadir))

%% Loading stimulus conditions

patterns   = {'rforidrft','rfsfdrft','posdisparitydrft','disparitydrft','cinterocdrft','coneinterocdrft','conedrft', ...
                'colorflicker','bwflicker','rfori','rfsize','cinteroc','color','rfsf','mcosinteroc','dotmapping'}; 

for p = 1:length(patterns)

   pattern      = patterns{p}; 
  
   if any(strfind(BRdatafile,pattern))
       startlog = strfind(BRdatafile,pattern); 
       if ~isequal(BRdatafile(startlog:end-3),pattern),continue
       else
       match    = patterns{p}; 
       end
   end
   
end

if isequal(match,'dotmapping')
ext  = '.gDotsXY_di';
else
ext  = ['.g' upper(match) 'Grating_di']; 
end

if contains(ext,'DRFT')
      grating     = readgDRFTGrating([filename ext]); % from nbanalysis 
elseif contains(ext,'Dots')
      grating     = readgDotsXY([filename ext]);
else
      grating     = readgGrating([filename ext]);
end

%% Load Event Times/Codes

NEV             = openNEV([filename '.nev'],'noread','overwrite');
EventCodes      = NEV.Data.SerialDigitalIO.UnparsedData - 128;        % we don't know why we have to subtract 128 but we do
EventSamples    = NEV.Data.SerialDigitalIO.TimeStamp;                 % in samples 
EventTimes      = floor(NEV.Data.SerialDigitalIO.TimeStampSec.*1000); % convert to ms 
[pEvC, pEvT]    = parsEventCodesML(EventCodes,EventSamples);          % sorts codes, samps or times into trials

% So far all of these data are from EVERY trial, including trials where
% animal breaks fixation. Lets get rid of those and make a new structure
% with the grating info and the stimulus onsets 

STIM            = sortStimandTimeData(grating,pEvC,pEvT,'stim'); % this is in nbanalysis. definitely double check it before you use it. 

%% STEP THREE: LOAD NEURAL DATA

%% LOAD LFP with NS2 file 
clear ext 

ext = 'ns2';
el  = 'eD';       % bank 
lfp = getLFP(filename,ext,el,'ascending');        % if you use this function you need to know the bank + sort direction 
                                                  % this function can be found in MLAnalysisOnline or nbanalysis
                                                  % if your sort direction input is correct, the data come out in order
                                                  % from top--> bottom
                                                  
%% quick way to get  bank info + sort direction if you don't have it in an Excel sheet or notes: 

% Read in NS Header

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

%% LOAD LFP and analog MUA with NS6 file ( you can follow these steps to load with ns2 as well but beware of sampling freq.)
%  let's break the whole thing down. 

clear ext NS_header banks neural 
% Read in NS Header
ext          = 'ns6'; 
NS_Header    = openNSx(strcat(filename,'.',ext),'noread');

% get basic info about recorded data
neural       = strcmp({NS_Header.ElectrodesInfo.ConnectorBank},el(2)); % logicals where contact bank name matches electrode of interest
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
    
    if neural(e) == 1    % why? because neural is a vector of logicals, 1 = contacts we want

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

%% Warning! THESE DATA ARE IN THE SAME ORDER AS THE BR PINS, NOT THE ORDER OF THE PROBE
%  sort data from top of electrode to bottom.
%  PLEASE LOOK AT THIS

sortdirection = 'ascending'; 

% sort electrode contacts in ascending order:
idx    = nan(length(NeuralLabels),1); 
for ch = 1:length(NeuralLabels)
    chname  = strcat(sprintf('%s',el),sprintf('%02d',ch));
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

STIM.LFP  = trigData(LFP,floor(STIM.onsets./30),-pre,post); % this function is MLAnalysisOnline or nbanalysis. pre variable is in absolute units 
STIM.CSD  = trigData(CSD,floor(STIM.onsets./30),-pre,post); 
STIM.aMUA = trigData(MUA,floor(STIM.onsets./30),-pre,post); 

%% LOAD discrete MUA with ppnev file
% to get the ppNEV files [post-processed NEV] use the offlineBRAutoSort
% directory under https://github.com/maierav/KiloSortUtils
% input is threshold (std) to apply to envelope to extract spikes

% let me show you where it is...

clear Fs 

load([filename '.ppNEV'],'-MAT','ppNEV')

Fs    = double(ppNEV.MetaTags.SampleRes);

for iii = 1:length(sortedLabels)
    clear elabel 
    elabel               = sortedLabels{iii}; 
    eidx                 = find(cell2mat(cellfun(@(x) contains(x',elabel),{ppNEV.ElectrodesInfo.ElectrodeLabel},'UniformOutput',0)));
    I                    = ppNEV.Data.Spikes.Electrode == eidx;
    SPK                  = double(ppNEV.Data.Spikes.TimeStamp(I)); % these are the spike times (in samples)
    sdf                  = spk2sdf(SPK,Fs); % this convolves the spikes. you need *jnm_kernel* to use the poisson dist 
   
    STIM.dMUA(:,iii,:)     = squeeze(trigData(sdf',floor(STIM.onsets./30),-pre,post));  % convolved ppNEV cmua sorted into trials 
    STIM.dMUA_bin(:,iii,:) = squeeze(trigBinaryData(SPK,pre,post,floor(STIM.onsets./30))); % binary ppNEV cmua sorted into trials 
end


%% LOOK at kilosorted data 
% use https://github.com/maierav/KiloSortUtils
% use phy to view results! 
% let me show you what it looks like on the MacPro...


%% STEP FOUR: plot the data
% LFP
basetp = STIM.LFP(25:75,:,:); %retrieves the reference range to which the data will be baseline corrected
mean_basetp = mean(basetp, 1); %calculates the mean of the reference range
bcd = STIM.LFP(1:351,:,:)- mean_basetp; %performs calculation, labeled bcd
mean_nLFP = mean(bcd, 3); %averages the contact samples by trial (the third dimension of the structure)
refwin = (pre:post); %reference window

h = figure('Position', [0,0,240,291*2]); %the dimensions of the figure
[ha, pos] = tight_subplot(24,1,[0.005 .03],[.1 .15],[.2 .2]); %channels, columns, [spacing], [top and bottom margin], [left and right margin]
for c = 1:24
    
    axes(ha(c)); % ha is a variable that gets the axis of each subplot, for each (i)1-24, get axis
    
    area(refwin,mean_nLFP(:,c)); %specify x range (refwin, mean_LFP(1 thru 351 samples, channel i)
    areafill = area(refwin,mean_nLFP(:,c));
    areafill.FaceColor = [0.75 0.75 0.75];
    xlim([-50 300]); %set the limitation of the x axis (sometimes the plot will plot more than you need, creating a gap)
    ylim([-320 100]);
    set(h,'position',get(h,'position').*[1 1 1 1]); %reduce the size of the individual plots
    line(xlim(), [0,0],'LineStyle','-.','LineWidth', 0.1, 'Color', 'k'); %create black horizontal line on the origin for each plot
    xline(0,'-.b');
    yticks(c);
    grid on
    
    if c < 24
    set(gca, 'XTick', '') %removing the units on the x axis (if i < 24) 
    %set(gca, 'YTick', '') %removing the units on the x axis (if i < 24)
    f = gca; f.XAxis.Visible = 'off'; %remove the x axis
    end
    
    if c == 1
        title({'\fontsize{16}',BRdatafile,'LFP'})
    end
   
    if c == 12
        ylabel('\fontsize{15}Contacts in order of depth');
    end
end

set(ha(1:23), 'XTickLabel',''); %remove labels from first 23 plots
set(ha(1:24), 'box', 'off');

xlabel('\fontsize{15}time (ms)');

%% aMUA
basetp = STIM.aMUA(25:75,:,:); %retrieves the reference range to which the data will be baseline corrected
mean_basetp = mean(basetp, 1); %calculates the mean of the reference range
bcd = STIM.aMUA(1:351,:,:)- mean_basetp; %performs calculation, labeled bcd
mean_naMUA = mean(bcd, 3); %averages the contact samples by trial (the third dimension of the structure)
refwin = (pre:post); %reference window

h = figure('Position', [0,0,240,291*2]); %the dimensions of the figure
[ha, pos] = tight_subplot(24,1,[0.005 .03],[.1 .15],[.2 .2]); %channels, columns, [spacing], [top and bottom margin], [left and right margin]
for c1 = 1:24
    
    axes(ha(c1)); % ha is a variable that gets the axis of each subplot, for each (c1)1-24, get axis
    
    area(refwin,mean_naMUA(:,c1)); 
    areafill = area(refwin,mean_naMUA(:,c1));
    areafill.FaceColor = [0.75 0.75 0.75];
    xlim([-50 300]); %set the limitation of the x axis (sometimes the plot will plot more than you need, creating a gap)
    ylim([-1 1]);
    set(h,'position',get(h,'position').*[1 1 1 1]); %reduce the size of the individual plots
    line(xlim(), [0,0],'LineStyle','-.','LineWidth', 0.1, 'Color', 'k'); %create black horizontal line on the origin for each plot
    xline(0,'-.b');
    %yticks(c1);
    ylabel(c1);
    grid on
    
    if c1 < 24
    set(gca, 'XTick', '') %removing the units on the x axis (if i < 24) 
    set(gca, 'YTick', '') %removing the units on the y axis (if i < 24)
    f = gca; f.XAxis.Visible = 'off'; %remove the x axis
    end
    
    if c1 == 1
        title({'\fontsize{16}',BRdatafile,'LFP'})
    end
   
    if c1 == 12
        ylabel({'\fontsize{13}Contacts in order of depth','\fontsize{9}12'})
    end
end

set(ha(1:23), 'XTickLabel',''); %remove labels from first 23 plots
set(ha(1:24), 'YTickLabel','');
set(ha(1:24), 'box', 'off');

xlabel('\fontsize{15}time (ms)');
%% Really useful Matlab central functions for plotting data to make everything easier and look better: 


% EXPORT FIG -- perfect for making pdfs when you are checking data!


% COLORS
% https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
% see also: http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3


% BAR PLOTS
% https://www.mathworks.com/matlabcentral/fileexchange/57499-superbarsuperbar 

