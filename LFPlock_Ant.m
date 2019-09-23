%% Written by Antonis Asiminas and Felicity Inkpen, SIDB July 2019.

% This script makes use of functions from FMAT, chronux, fieldtrip and PAC toolbox as well as some functions from Kitanishi et al 2015 Neuron 
% This was written for Asiminas et al 2019 eLife
% The inputs are: 
% mypos position data
% .set file from Axona system
% .egf files recorded in Axona system which contain 4800Hz continuous LFP data, 
% an xl file containing cluster identities from previously analysed data,
% .clu and .fet files from Klusters were I get the cluster spike times.
% These are from previsouly combined sessions in order to spike sort the same clusters accross them.  

% The clusters used in this analysis have been chosen because they are
% putative pyramidal cells nicely isolated from noise. They are the same
% clusters that have been use in SpikeTrainAnt

%The way I run this is pretty silly. I run it separatelly for every
%rat/session and then get the outputs, average across genotypes and run
%stats in SPSS or R and create graphs in Graphpad/Inkscape


% The script is separated into four sections:
% In the fist section the position data are analysed and periods of below threshold mobility are
% identified
% In the second section the LFP data are imported and analysed. Basic
% spectrograms and the power of oscillatory bands of interest are
% calculated. Phase-amplitude coherences are computed and their recurrence
% during the duration of the experiment is explored by correlations
% In the third section the spike timestamps are imported and are corrected
% if they're coming from previously combined sessions (session number >1)
% In the fourth section the LFP and spike timestamps are coming together in
% order to calculate phase locking of spikes to different oscillatory
% bands (theta, low gamma and high gamma). 
% The last section is about taking everything and creating plots and output tables for further analysis.

%% Inputs and parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; %closes open figures and files
clear all; %clears workspace from previous runs of the script
disp('Starting analysis...');

sessions = input('How many sessions(number without quotes i.e 6)?');
SesDurMin = input('How many minutes was each session?');
experiment = input('Base name of sessions(in single quotes i.e ''20160516'')?');

tic;
elecs = [1 2 3 4 5 6 7 8]; %I don't really need this vector but it makes my life easier for the loops when I read clusters later on

%Limits that wil be used in fft later in the script
%These are based on Kitanishi et al 2015 Neuron but also very close to
%Dvorak et al 2018 PLoS Biol. I call it fast gamma but technically is
%medium gamma
% I could have done that with a struct but it's difficult to index if I
% need to
% Theta bandpass
Lims{1,1} = 5; %low stop for theta
Lims{1,2} = 6; %low pass for theta
Lims{2,1} = 11; %high pass for theta
Lims{2,2} = 12; %high stop for theta
% Slow gamma bandpass
Lims{1,3} = 18; %low stop for slow gamma
Lims{1,4} = 20; %low pass for slow gamma
Lims{2,3} = 48; %high pass for slow gamma
Lims{2,4} = 50; %high stop for slow gamma
% Fast gamma bandpass
Lims{1,5} = 58; %low stop for fast gamma
Lims{1,6} = 60; %low pass for fast gamma
Lims{2,5} = 88; %high pass for fast gamma
Lims{2,6} = 90; %high stop for fast gamma

% Gamma event window length [sec]
% Spikes which happend in (gamma power peak)+-(wLength/2) will be used.
% Default, 0.400 sec according to Laura et al., 2009.
wLength = 0.400;
pixelratio=330; %that's pixels per meter. This has to be adjusted for the open ephys system recordings.

startup; %this call a function at the bottom of the script for figure uniformity

%-------------------------------------%%%%%%%%%%%%%%%%%
%% Section One %%%% Loading position data that will be used for excluding spikes and LFP during imobility.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Reading positions ...');

pos_all=load('merge.mypos'); %I'm loading the merged position file. This is smoothed compared to the individual session position files.

%Extract series of coordinates
pos_x_all = pos_all(:,1);
pos_y_all = pos_all(:, 2);
pos_time_all = pos_all(:,4)*100000; %that's just to be able to use the Tseg file to split
disp('...done');

T_segs =load('merge.tseg'); %use the tseg file to see the size of the segments and the number of them
T_segs=T_segs(:,2);
T_segs = [T_segs(1:end);T_segs(end)+(T_segs(end)-T_segs(end-1))]*100000; %that is to correct with the spikes later

%this separates the x-t-t data into sessions
disp('Analysing movement...');
Pos_X_all_ses = zeros(SesDurMin*60*50, sessions); %that's based on the assumption that all session should be the same length. Probably not the most elegant way to write that. If they were not the same I should use a cell array anyway
Pos_Y_all_ses = zeros(SesDurMin*60*50, sessions);
Pos_Time_all_ses = zeros(SesDurMin*60*50, sessions); %I need 600sec*50Hz = 30000 samples for each of the sessions. Probably not the most elegant way to write that.

for tp=1:length(T_segs)-1 
    vals=[]; %clean vals after every iteration
    ind_start=find(pos_time_all==T_segs(tp));
    ind_dif=find(pos_time_all==T_segs(tp+1))-find(pos_time_all==T_segs(tp));
    if tp==6
        ind_dif=length(pos_time_all)-find(pos_time_all==T_segs(tp));
    end
    
    if ind_dif>30000  %usually a 10min session in axona is about 30049 datapoints
        n=1;
        while n<SesDurMin*60*50+1    %here I take only the datapointsthat correspond to 50Hz*Session in seconds
        vals(n,1) = pos_time_all(ind_start);
        ind_start=ind_start+1;
        n=n+1;
        end 
        
    Pos_Time_all_ses(:,tp)=vals-T_segs(tp); %If everything is right then all columns of Pos_Time_all_ses should be the same. I also realised that the system has a lag of about 1sec (0.96sec). That means that for that period the EEG is 0 and nothing is happening. I want correct that by taking only the session time. So I need 600sec*50Hz = 30000 samples for each of these   
    [~,inds]=intersect(pos_time_all, vals); %find the indeces of each separate session smpaling times in the over alltime vector
    Pos_X_all_ses(:,tp)=pos_x_all(inds);  
    Pos_Y_all_ses(:,tp)=pos_y_all(inds);
    
    else %sometimes a session will be a few datapoints short (usually 29999 for the 10min sessions I have analysed before). 
        disp(['Session ', num2str(tp),'has only ',num2str(ind_dif), 'position data points'])
        
        n=1;
        while n<ind_dif %with this I take all the available time datapoints ...
        vals(n,1) = pos_time_all(ind_start);
        ind_start=ind_start+1;
        n=n+1;
        end
        
        [~,inds]=intersect(pos_time_all, vals);%... and the corresponding X and Y coordinates...
        X_help= pos_x_all(inds);
        Y_help= pos_y_all(inds);
        
        
        while length(vals)<SesDurMin*60*50 %...and I complete the remaining few datapoints with fake X-Y coordinates
            vals(end+1)=vals(end)+2000;
            X_help(end+1)=X_help(end);
            Y_help(end+1)=Y_help(end);
        end
     
     Pos_Time_all_ses(:,tp)=vals-T_segs(tp);
     Pos_X_all_ses(:,tp)= X_help;  
     Pos_Y_all_ses(:,tp)= Y_help;  
       
    end
   
end

%% I'm currently only using that for Fliss' bit not for the rest of my analysis because the sliding windows make the indexing between electrophysiological data and speed difficult.
%I know this is counterintuitive but I will try to tidy up the script on a later phase
% Checking speed of movement for all sessions. To do that I'm using 500ms
% bins that slide (step 40ms).


slwindow=25; %this is based on the sampling rate every 0.02s (20ms) 50Hz.Therefore 25 samples are 500ms
step=2; %this makes the sliding window 0.04s (40ms). That's a 12.5% overlap. Trimper used 50ms and Arbab 25ms

Velocity_cmpersec_all=zeros(floor((length(Pos_Time_all_ses)-slwindow)/2),sessions); %this is a matrix containing the speed in cm/sec for all the sliding 500ms windows (rows) for all sessions (columns). I assume that all sessions have the same length so 
Tsegs_stationary=cell(1,sessions); %this is a cell array containing the 500ms time segments the animal was below 3cm/s for all sessions. This will be use to exclude periods for our LFP analysis and exclude spikes in the spike phase locking bit of the script
time_lim_vec=[]; 
for h=1:sessions
    time_lim_vec=zeros(floor((length(Pos_Time_all_ses)-slwindow)/2),2); %this contains all the start and end of all 500ms bins for which I'll check the speed
    e=1;
    for g=1:step:length(Pos_Time_all_ses(:,h))-slwindow %the last time point is cropped because of the sliding window but that not a problem.
        time_lim_vec(e,1)=Pos_Time_all_ses(g,h);
        time_lim_vec(e,2)=Pos_Time_all_ses(g+slwindow,h); 
        e=e+1;
    end 

    for g=1:length(time_lim_vec) %for all the 500ms segments (with 460ms overlap)
        start_bin=find(Pos_Time_all_ses(:,h)==time_lim_vec(g,1)); %find the index of the start of this segment on the Sessions cell
        end_bin=find(Pos_Time_all_ses(:,h)==time_lim_vec(g,2)); %find the index of the end of this segment on the Sessions cell
        dist_bin=((sqrt((Pos_X_all_ses(end_bin,h)-Pos_X_all_ses(start_bin,h)).^2 + (Pos_Y_all_ses(end_bin,h)-Pos_Y_all_ses(start_bin,h)).^2))*100)/pixelratio; %calculate distance (Pythagoras theorem, Go Samos!) traveled and convert pixels to cm based on ratio. The 100 is because the pixelratio is in pixels per meter
        vel_bin=(dist_bin/(time_lim_vec(g,2)-time_lim_vec(g,1)))*100000; %calculate the velocity for that segment. I convert time back to seconds 
        Velocity_cmpersec_all(g,h)=vel_bin; %add it to overview matrix
    end
   
   notmove_start=time_lim_vec(Velocity_cmpersec_all(:,h)<3 | Velocity_cmpersec_all(:,h)>80, 1); %the start of the window with velocity less than 3cm/s OR unreasonably fast
   notmove_end=time_lim_vec(Velocity_cmpersec_all(:,h)<3 | Velocity_cmpersec_all(:,h)>80, 2); %the end of the window with velocity less than 3cm/s OR unreasonably fast
   Tsegs_stationary{h}(:,1)= notmove_start;
   Tsegs_stationary{h}(:,2)= notmove_end;
end


%--------------------------------%%%%%%%%%%%%%%%%%%
%% Section Two %%%%%%%%%%%%%%%%% Reading EGF file and doing some basic analyses of the LFP data. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This section will provide a number of outputs suitable for plotting but
% also some variables that will be used together with the spiketiming data
% (Section Three)for the phase locking analysis.
% Some of fuctions that will be used here are from the PAC toolbox (Onslow,
% May 2010)


%Reading EGF files and their sampling frequencies from the folder and
%adding them in cell arrays
EGF_all = cell(1,sessions); %a cell array with the EEG for all sessions
Fs_all = zeros(1,sessions); %a vector with the sampling rate of all sessions' EEG
for ses=1:sessions
    egffile = [experiment, num2str(ses),'.egf'];
    disp(['Reading eeg data from session ', num2str(ses)])
    [EGF,Fs] = readEGF(egffile); %function at the bottom of the script
    if length(EGF)< Fs*60*SesDurMin+ 4800   %The eeg files are supposed to have 1 more second worth of samples that are 0. Fliss' power spectral density analysis section needs that. If it's not there I'll just add it.
        
        while length(EGF)-(Fs*60*SesDurMin + 4800)<0 %for my analysis later (line 388) I only take the samples for Fs*60*SesDurMin (for a 10min session in 4.8kHz sampling rate it's 2880000)
            EGF(end+1)=0;     
        end
    end
      
    EGF_all{ses}=EGF;
    Fs_all(ses)=Fs;
end


%Depending on the way LFP was recorded in the axona system, it might need
%to be inverted. For that we use this EegMode function to see what is the
%mode and also see the gains in that channels so we can calculate the mV.
%We get that info from the .set file
%a cell array with the mode of all sessions' EEG. This might be useful for checking later I used cell arrays because for some rats I might have recorded eeg from more than one channel
eegmode_all = zeros(1,sessions); 
Gain_all = zeros(1,sessions); %a cell array with the gains sessions' EEG. This will be used to get ?V later.
for ses=1:sessions
    setfile = [experiment, num2str(ses),'.set'];
    [eegmode, Gain] = EegMode(setfile);
    if eegmode == 3
                EGF_all{ses} = -EGF_all{ses};
                disp(['Recorded in -Ref mode! EEG sign was corrected in session', num2str(ses)])
    end
    eegmode_all(ses)=eegmode;
    Gain_all(ses)=Gain;
end


% Convert EEG values to micro Volt (uV)
EGF_all_uV=cell(1,sessions); %a cell array with the EEG for all sessions in micro Volt values
for ses=1:sessions
EGF_all_uV{ses} = EGF_all{ses} ./ 32768 .* 1500000 ./ Gain_all(ses); % these are stored at 16-bit resolution, so the values range from -32768 to 32767
end







% Raw text from Jim outlining the process above:

%     The values in the files are in bits.  That is, they have not been adjusted
%     for the gain factor.  To convert those numbers to microvolts, you need to
%     do three things:
%
%     1. Look up the gain for the channel in the .set file for the trial.  The
%     relevant values are called gain_ch_X, where X goes from 0 to 127 for
%     channels 1 to 128.  So, if you want to convert the samples for channel 1
%     to volts, look up gain_ch_0.
%
%     2. Next, look up the value of ADC_fullscale_mv in the .set file.  It is
%     almost always 1500 for dacqUSB systems.
%
%     3. The data values for .egf files range from -32768 to 32767 (16 bits).
%
%     To convert from sample value to microvolts, the factor is:
%
%     value_in_uV = sample_value / 32768 * ADC_fullscale_mv * 1000 / gain
%                                    (ADC_fullscale_mv 1500)
%     where "sample_value" is the signed-integer value you read from the file,
%     "ADC_fullscale_mv" is the value you found at step 2 above, and "gain" is
%     the relevant gain factor you found at step 1.
%
%     For spike (.N) or .eeg data, these are stored at 8-bit resolution, so the
%     values range from -128 to 127, and the conversion factor is:
%
%     value_in_uV = sample_value / 128 * ADC_fullscale_mv * 1000 / gain



%this will give me the duration of each session. If everything is ok the
%durations should be the same (601sec) for every session. This will feed
%into the main function that brings LFP and spike timing together.
Ses_Dur_all=zeros(1,sessions);
for ses=1:sessions
    setfile = [experiment, num2str(ses),'.set'];
    Ses_Dur = SessionDuration(setfile);
    Ses_Dur_all(ses)=Ses_Dur;
end

%% --Calculation of spectrograms from Fliss. This uses some previously calculated variables and calculates power density in a number of different ways. At this point only a subset of these is exported
%Create a vector of values to normalise the spectrogram data to in the decibel conversion
[baseline_normalisation]=baseline_norm(Velocity_cmpersec_all, Fs, Fs_all, sessions, EGF_all_uV);

%The function used below is from the FMAT toolbox and uses stuff from the chronux
%toolbox (I know it's annoying to track down shit like that).It has two more possible outputs but I don't know what to do with them. It plots a
%spectrogram but it's across time 600sec (10min session). Colour indicates
%power but not easy to see/compare
Spectro_all =cell(1, sessions);
db_all =cell(1, sessions);
lgtrfm_all=cell(1, sessions);
for ses=1:sessions
    t = (0:length(EGF_all_uV{ses})-1)'/Fs_all(ses); %I need to create time-voltage pairs for this function.
  %  [spectrogram,t_spect,f_spect] = MTSpectrogram([t EGF_all_uV{ses}],'show','on','frequency',Fs_all(ses), 'range',[0 100]); %I used the lowest and the highest frequency I'm interested in this.
    [spectrogram, logTransformed, dbconverted,t_spect,f_spect] = MTSpectrogram([t EGF_all_uV{ses}],baseline_normalisation,'show','on','frequency',Fs_all(ses), 'range',[0 100]); %I used the lowest and the highest frequency I'm interested in this.

    Spectro_all{ses}=spectrogram;
    db_all{ses}=dbconverted;
    lgtrfm_all{ses}=logTransformed;

figure
subplot(1,2,1)
plot(f_spect,spectrogram(:,1));
yyaxis right
plot(f_spect,logTransformed(:,1));
legend('LFP data', 'Log transformed data')
legend boxoff
title('Power-frequency spectra, example trace')

subplot(1,2,2)
plot(f_spect,spectrogram(:,1));
yyaxis right
plot(f_spect,dbconverted(:,1)); 
legend('LFP data', 'dB corrected data')
legend boxoff
title('Power-frequency spectra, example trace')



figure
subplot(1,2,1)
plot(f_spect,mean(spectrogram,2));
yyaxis right
plot(f_spect,mean(logTransformed,2));
legend('LFP data', 'Log transformed data')
legend boxoff
title('Power-frequency spectra, mean trace')

subplot(1,2,2)
plot(f_spect,mean(spectrogram,2));
yyaxis right
plot(f_spect,mean(dbconverted,2));
legend('LFP data', 'dB corrected data')
legend boxoff
title('Power-frequency spectra, mean trace')
end

%Create minute-by-minute mean frequency spectra
minute_by_minute_spectra(sessions, f_spect, Spectro_all)

%observe the mean power spectra per session with standard deviation
boundedline_powerspectra(sessions, Spectro_all,f_spect);
suptitle (['Experiment ' experiment ', power spectra'])
boundedline_powerspectra(sessions, lgtrfm_all,f_spect);
suptitle (['Experiment ' experiment ', log transformed power spectra'])
boundedline_powerspectra(sessions, db_all,f_spect);
suptitle (['Experiment ' experiment ', decibel converted power spectra'])

%create a normalised-to-the-data version of the power spectra, normalised
%to the mean value between 4 and 100 Hz
for i=1:length(Spectro_all)
Spectro_all_norm{1,i}=mean(Spectro_all{1,i},2);
end
normalisation_factor=cell2mat(Spectro_all_norm);
normalisation_factor=max(max(normalisation_factor(28:end,:)));
for ses=1:sessions
Spectro_all_norm{1,ses}=Spectro_all_norm{1,ses}./normalisation_factor;
end
figure
boundedline_normalised(sessions, Spectro_all_norm,Spectro_all,f_spect, normalisation_factor, experiment);

% Averaged spectrogram output

% function to call within the script. Should occur after the boundedline_normalised function is called. 

[spectrogram_output]=average_spectrogram(Spectro_all,sessions, lgtrfm_all,db_all,Spectro_all_norm,normalisation_factor);


% calling the function, to be included at the end of the section I have worked on 

[frequency_bands, theta_spect, slow_gamma_spect, fast_gamma_spect]=frequency_bands_spectra(f_spect, Spectro_all, sessions, SesDurMin);


close all %this is to close all generated figures from Fliss' section
%% This is to reshape the long EEG signal of every session to 500ms non-overlapping epochs. I also take only the first 600sec for each session and not the lag 0.93s at the end
% From these segments I can export some examples. I just need to run the fft bandpass filters on a good example section(s)

EGF_500ms_all=cell(1, sessions);
for ses=1:sessions
EGF_500ms=reshape(EGF_all_uV{ses}(1:Fs*60*SesDurMin), [Fs*0.5 60*SesDurMin/0.5]); 
EGF_500ms_all{ses}=EGF_500ms;
end

% ** See line 384. I find the cells that have the highest Theta power and
% the highest slow gamma and high gamma power and choose examples
%% Thresholding all these 500ms non-overlapping segments for theta power 
disp('Computing how LFP power is modulated by animal speed')
Time_ses=(0:length(EGF_all_uV{ses}(1:Fs*60*SesDurMin))-1)'/Fs +1/Fs; %that also corrects for the time of the first sampling with should be 1/Fs
Time_500ms_all=reshape(Time_ses,[Fs*0.5 1200]); %this will be used to index time segments. =< from first row and >from last row of each segment (column)
%Some parameters I need for pwelch function. 
%To be honest this is a bit stupid since Fliss is working on spectrogram and band powers but I need these to threshold based on theta power and also check how speed correlates with the power of different bands.
%In a future, cleaner version, I should try to calculate things once and use them throughout.
freqRange=Lims{1,1}:0.5:Lims{2,6}; %range of frequencies in 0.5 Hz steps. That's within the frequencies of interest 
segmentLength=[]; noverlap=[]; %default Hamming (8 windows with 50% overlap)

PWR_bands_500_all=cell(1,sessions); %this has the theta, low gamma and high gamma band power for every 500ms non-overlapping epoch in each session. I can plot them against speed or sort based on theta power

for ses=1:sessions
    [Pxx,F] = pwelch(EGF_500ms_all{ses},segmentLength,noverlap,freqRange,Fs); %pwelch calculates PSD independently for each column (500ms bin) and stores in corresponding column of Pxx
    
    for epoch500=1:size(Pxx,2) %for every epoch of 500ms caclulate the power of theta and put it in a vector
        thetapower=bandpower(Pxx(:,epoch500),F,[Lims{1,2} Lims{2,1}],'psd');
        lowgammapower=bandpower(Pxx(:,epoch500),F,[Lims{1,4} Lims{2,3}],'psd');
        highgammapower=bandpower(Pxx(:,epoch500),F,[Lims{1,6} Lims{2,5}],'psd');
       
        PWR_bands_500_all{ses}(1,epoch500)=thetapower;
        PWR_bands_500_all{ses}(2,epoch500)=lowgammapower;
        PWR_bands_500_all{ses}(3,epoch500)=highgammapower;
        
    end
end

%Normalise the bands to the average power of the whole session for that band. 
PWR_bands_500_all_norm=cell(1,sessions); 

for ses=1:sessions
    for band=1:length(Lims)/2
        PWR_bands_500_all_norm{ses}(band,:)=PWR_bands_500_all{ses}(band,:)/mean(PWR_bands_500_all{ses}(band,:),'all');
    end
end

%Also generating the LG/HG power; this is to explore data in similar way to Dvorak et al 2018 PLoS Biol
%I NEED TO EXPORT THIS MATRIX AND THE MATRIX Velocity_cmpersec_500 SO I CAN
%PLOT THIS POWER RATIO AGAINST VELOCITY
PWR_LG2HG_500_ratio=NaN(sessions,size(PWR_bands_500_all{1}, 2)); %this has the LG/HG ratio of every session(rows) for every 500ms epoch (columns)
for ses=1:sessions

    PWR_LG2HG_500_ratio(ses,:)=PWR_bands_500_all_norm{ses}(2,:)./PWR_bands_500_all_norm{ses}(3,:);
end


% This is to find some nice example traces for figures. I NEED TO PLOT
% THESE EXAMPLES AS SVG AND THEN USE AN EXAMPLE FROM EACH SESSION IN FIGURE
% 3 OF THE PAPER
% I'm finding the epochs when I have nice theta and gamma power
StrongLFP_bandpass=cell(1+length(Lims)/2,sessions); %this will have 4 example traces for each session from that rat. Top is raw, second is theta, third is slow gamma and fourth is fast gamma
for ses=1:sessions
    
    [Tq1, Tq2, Tq3] = findiquartile(PWR_bands_500_all{ses}(1,:));
    [Sq1, Sq2, Sq3] = findiquartile(PWR_bands_500_all{ses}(2,:));
    [Fq1, Fq2, Fq3] = findiquartile(PWR_bands_500_all{ses}(3,:));
    
    Ttop= PWR_bands_500_all{ses}(1,:)>=Tq1; %depending what the findings from the band power are, I can change that a get a medium power theta or gamma to illustrate a point
    Stop= PWR_bands_500_all{ses}(2,:)>=Sq1;
    Ftop= PWR_bands_500_all{ses}(3,:)>=Fq1;
    
    Alltop = Ttop & Stop & Ftop; %That is a logical vector. Ones indicate the 500ms epochs where theta, slow gamma and high gamma are strong.
    
    StrongLFP=EGF_500ms_all{ses}(:,Alltop);
    indx=2;
    for i=1:2:length(Lims)-1
        %I'm taking the middle segment ceil(end/2) from all the segments
        %that meet the criterion in line 451
        StrongLFP_bandpass{1,ses}=StrongLFP(:,ceil(end/2)); %getting the raw LFP signal
        StrongLFP_bandpass{indx,ses} = (fftbandpass(StrongLFP(:,ceil(end/2)),Fs, Lims{1,i}, Lims{1,i+1}, Lims{2,i}, Lims{2,i+1}))'; % Theta, low-gamma and high-gamma bands         
        indx=indx+1;
    end
    
end



% Calculate speed for all the epochs. I can export that for the movement analysis to show that genotypes did not different in speed of movement.
% Then remove the bins that correspond to EEG outliers so I can plot speed against band power
Velocity_cmpersec_500=NaN(sessions,size(EGF_500ms,2)); %EGF_500ms is one of the matrices stored in the cell array EGF_500ms_all. All of them have the same size with number of columns being 600s/0.5s =1200
for ses=1:sessions
    for epoch=1:size(Time_500ms_all, 2)
        movX_all=Pos_X_all_ses(Pos_Time_all_ses(:,ses)>=Time_500ms_all(1,epoch)*100000 & Pos_Time_all_ses(:,ses)<=Time_500ms_all(end,epoch)*100000,ses); %getting ll the  X,Y points during the 500ms epoch of the iteration
        movY_all=Pos_Y_all_ses(Pos_Time_all_ses(:,ses)>=Time_500ms_all(1,epoch)*100000 & Pos_Time_all_ses(:,ses)<=Time_500ms_all(end,epoch)*100000,ses); 
        startX=movX_all(1); %getting the starts and ends of for X,Y movement for each 500ms epoch 
        endX=movX_all(end);
        startY=movY_all(1);
        endY=movY_all(end);
        Velocity_cmpersec_500(ses,epoch)= (((sqrt((endX-startX)^2 + (endY-startY)^2))*100)/pixelratio)/(Time_500ms_all(end,epoch)-Time_500ms_all(1,epoch)); %that the speed for that epoch
    end
end


% Removing outliers from all bands in PWR_bands_500_all and PWR_LG2HG_500_ratio
%I NEED TO EXPORT THE MATRICES AND VECTORS FROM BOTH
%PWR_bands_500_nooutliers AND Vel_4_PWR_bands SO I CAN PLOT THE MODULATION
%OF LFP POWER BY ANIMALS VELOCITY.
PWR_bands_500_nooutliers=cell(length(Lims)/2, sessions); 
Vel_4_PWR_bands=cell(length(Lims)/2, sessions); %I have different velocity matrices for the different bands because I need to exclude outliers for every band and that could mean that not all bands will then have the same 500ms epochs included 
PWR_LG2HG_500_ratio_nonoutliers=cell(1,sessions);
Vel_4_LG2HG_500_ratio=cell(1,sessions);
for ses=1:sessions%for all sessions
    
    PWR_LG2HG_500_ratio_nonoutliers{ses}=PWR_LG2HG_500_ratio(ses,(~isoutlier(PWR_LG2HG_500_ratio(ses,:), 'median')));
    Vel_4_LG2HG_500_ratio{ses} = Velocity_cmpersec_500(ses,(~isoutlier(PWR_LG2HG_500_ratio(ses,:), 'median')));
    for band=1:length(Lims)/2 %and all frequency bands in PWR_bands_500_all       
        PWR_bands_500_nooutliers{band, ses}=PWR_bands_500_all_norm{ses}(band,(~isoutlier(PWR_bands_500_all_norm{ses}(band,:), 'median'))); %removing outlier for all three bands of interest
        Vel_4_PWR_bands{band,ses}= Velocity_cmpersec_500(ses,(~isoutlier(PWR_bands_500_all_norm{ses}(band,:), 'median'))); %removing speed for the same outliers for all three bands of niterest
    end
    
end

%% Here I calculate the PAC for different powers of Theta in all sessions. I won't do it in a minute by minute at this stage. I have commented out that further down.

%Calculating Phase-Amplitude coherence. The method used is the Modulation Index (Canolty et al 2006 Science)
Time_IQi_THPWR_all=cell(4,sessions); %this is to capture the start and end of the 500ms epochs that theta power is in the 1st, 2nd, 3rd or 4th interquartile (IQ) interval
PAC_HPCt2lg_all=cell(4,sessions); %these are to capture the PAC, PACmean and PACmax low gamma and fast gamma against theta for all sessions at all 4 interquartile (IQ) intervals of theta power
PAC_HPCt2hg_all=cell(4,sessions);
PAC_HPCt2lg_means=zeros(4,sessions);
PAC_HPCt2lg_max=zeros(4,sessions);
PAC_HPCt2hg_means=zeros(4,sessions);
PAC_HPCt2hg_max=zeros(4,sessions);
disp('Computing Phase-Amplitute Coupling (PAC) for all sessions and for four theta power interquartile (IQ) intervals. This will take some time...')
for ses=1:sessions
    
    [q1, q2, q3] = findiquartile(PWR_bands_500_all{ses}(1,:)); % see info at the end of the script. This is based on Tort et al 2009 PNAS (Fig 4 B more specifically)
    
    
    Time_1IQi_THPWR(:,1)=Time_500ms_all(1,(PWR_bands_500_all{ses}(1,:)<q1)); %with these few few lines I find the start and the end of epochs that have a theta power in the 1st, 2nd, 3rd or 4th interquartile (IQ) interval
    Time_1IQi_THPWR(:,2)=Time_500ms_all(end,(PWR_bands_500_all{ses}(1,:)<q1));
    
    %PAC between theta and low gamma for 1st interquartile (IQ) interval of ThetaPWR
    [PAC_HPCt2lg_1IQi,~, ~]=find_pac_shf(EGF_500ms_all{ses}(:,(PWR_bands_500_all{ses}(1,:)<q1)),Fs,'mi',EGF_500ms_all{ses}(:,(PWR_bands_500_all{ses}(1,:)<q1)),Lims{1,2}:0.5:Lims{2,1},Lims{1,4}:0.5:Lims{2,3},'n',0,7,200,0,0.05,'HPC-HPC','HPC','HPC');
    %PAC between theta and high gamma for 1st interquartile (IQ) interval of ThetaPWR
    [PAC_HPCt2hg_1IQi,~, ~]=find_pac_shf(EGF_500ms_all{ses}(:,(PWR_bands_500_all{ses}(1,:)<q1)),Fs,'mi',EGF_500ms_all{ses}(:,(PWR_bands_500_all{ses}(1,:)<q1)),Lims{1,2}:0.5:Lims{2,1},Lims{1,6}:0.5:Lims{2,5},'n',0,7,200,0,0.05,'HPC-HPC','HPC','HPC');
    
    
    Time_2IQi_THPWR(:,1)=Time_500ms_all(1,(PWR_bands_500_all{ses}(1,:)>=q1 & PWR_bands_500_all{ses}(1,:)<q2));
    Time_2IQi_THPWR(:,2)=Time_500ms_all(end,(PWR_bands_500_all{ses}(1,:)>=q1 & PWR_bands_500_all{ses}(1,:)<q2));
    
    %PAC between theta and low gamma for 2nd interquartile (IQ) interval of ThetaPWR
    [PAC_HPCt2lg_2IQi,~, ~]=find_pac_shf(EGF_500ms_all{ses}(:,(PWR_bands_500_all{ses}(1,:)>=q1 & PWR_bands_500_all{ses}(1,:)<q2)),Fs,'mi',EGF_500ms_all{ses}(:,(PWR_bands_500_all{ses}(1,:)>=q1 & PWR_bands_500_all{ses}(1,:)<q2)),Lims{1,2}:0.5:Lims{2,1},Lims{1,4}:0.5:Lims{2,3},'n',0,7,200,0,0.05,'HPC-HPC','HPC','HPC');
    %PAC between theta and high gamma for 2nd interquartile (IQ) interval of ThetaPWR
    [PAC_HPCt2hg_2IQi,~, ~]=find_pac_shf(EGF_500ms_all{ses}(:,(PWR_bands_500_all{ses}(1,:)>=q1 & PWR_bands_500_all{ses}(1,:)<q2)),Fs,'mi',EGF_500ms_all{ses}(:,(PWR_bands_500_all{ses}(1,:)>=q1 & PWR_bands_500_all{ses}(1,:)<q2)),Lims{1,2}:0.5:Lims{2,1},Lims{1,6}:0.5:Lims{2,5},'n',0,7,200,0,0.05,'HPC-HPC','HPC','HPC');
    
    
    Time_3IQi_THPWR(:,1)=Time_500ms_all(1,(PWR_bands_500_all{ses}(1,:)>=q2 & PWR_bands_500_all{ses}(1,:)<q3));
    Time_3IQi_THPWR(:,2)=Time_500ms_all(end,(PWR_bands_500_all{ses}(1,:)>=q2 & PWR_bands_500_all{ses}(1,:)<q3));
    
    %PAC between theta and low gamma for 3nd interquartile (IQ) interval of ThetaPWR
    [PAC_HPCt2lg_3IQi,~, ~]=find_pac_shf(EGF_500ms_all{ses}(:,(PWR_bands_500_all{ses}(1,:)>=q2 & PWR_bands_500_all{ses}(1,:)<q3)),Fs,'mi',EGF_500ms_all{ses}(:,(PWR_bands_500_all{ses}(1,:)>=q2 & PWR_bands_500_all{ses}(1,:)<q3)),Lims{1,2}:0.5:Lims{2,1},Lims{1,4}:0.5:Lims{2,3},'n',0,7,200,0,0.05,'HPC-HPC','HPC','HPC');
    %PAC between theta and high gamma for 3nd interquartile (IQ) interval of ThetaPWR
    [PAC_HPCt2hg_3IQi,~, ~]=find_pac_shf(EGF_500ms_all{ses}(:,(PWR_bands_500_all{ses}(1,:)>=q2 & PWR_bands_500_all{ses}(1,:)<q3)),Fs,'mi',EGF_500ms_all{ses}(:,(PWR_bands_500_all{ses}(1,:)>=q2 & PWR_bands_500_all{ses}(1,:)<q3)),Lims{1,2}:0.5:Lims{2,1},Lims{1,6}:0.5:Lims{2,5},'n',0,7,200,0,0.05,'HPC-HPC','HPC','HPC');
    
    
    Time_4IQi_THPWR(:,1)=Time_500ms_all(1,(PWR_bands_500_all{ses}(1,:)>=q3));
    Time_4IQi_THPWR(:,2)=Time_500ms_all(end,(PWR_bands_500_all{ses}(1,:)>=q3));
    
    %PAC between theta and low gamma for 4th interquartile (IQ) interval of ThetaPWR
    [PAC_HPCt2lg_4IQi,~, ~]=find_pac_shf(EGF_500ms_all{ses}(:,(PWR_bands_500_all{ses}(1,:)>=q3)),Fs,'mi',EGF_500ms_all{ses}(:,(PWR_bands_500_all{ses}(1,:)>=q3)),Lims{1,2}:0.5:Lims{2,1},Lims{1,4}:0.5:Lims{2,3},'n',0,7,200,0,0.05,'HPC-HPC','HPC','HPC');
    %PAC between theta and high gamma for 4th interquartile (IQ) interval of ThetaPWR
    [PAC_HPCt2hg_4IQi,~, ~]=find_pac_shf(EGF_500ms_all{ses}(:,(PWR_bands_500_all{ses}(1,:)>=q3)),Fs,'mi',EGF_500ms_all{ses}(:,(PWR_bands_500_all{ses}(1,:)>=q3)),Lims{1,2}:0.5:Lims{2,1},Lims{1,6}:0.5:Lims{2,5},'n',0,7,200,0,0.05,'HPC-HPC','HPC','HPC');
    
    %capture everything I need for Theta-Low gamma
    PAC_HPCt2lg_all{1,ses}=PAC_HPCt2lg_1IQi;
    PAC_HPCt2lg_all{2,ses}=PAC_HPCt2lg_2IQi;
    PAC_HPCt2lg_all{3,ses}=PAC_HPCt2lg_3IQi;
    PAC_HPCt2lg_all{4,ses}=PAC_HPCt2lg_4IQi;   
    
    PAC_HPCt2lg_means(1,ses)=mean(PAC_HPCt2lg_1IQi,'all');
    PAC_HPCt2lg_means(2,ses)=mean(PAC_HPCt2lg_2IQi,'all');
    PAC_HPCt2lg_means(3,ses)=mean(PAC_HPCt2lg_3IQi,'all');
    PAC_HPCt2lg_means(4,ses)=mean(PAC_HPCt2lg_4IQi,'all');
    
    PAC_HPCt2lg_max(1,ses)=max(PAC_HPCt2lg_1IQi,[],'all');
    PAC_HPCt2lg_max(2,ses)=max(PAC_HPCt2lg_2IQi,[],'all');
    PAC_HPCt2lg_max(3,ses)=max(PAC_HPCt2lg_3IQi,[],'all');
    PAC_HPCt2lg_max(4,ses)=max(PAC_HPCt2lg_4IQi,[],'all');
    
    %--------------- for Theta-High gamma
    
    PAC_HPCt2hg_all{1,ses}=PAC_HPCt2hg_1IQi;
    PAC_HPCt2hg_all{2,ses}=PAC_HPCt2hg_2IQi;
    PAC_HPCt2hg_all{3,ses}=PAC_HPCt2hg_3IQi;
    PAC_HPCt2hg_all{4,ses}=PAC_HPCt2hg_4IQi;
      
    PAC_HPCt2hg_means(1,ses)=mean(PAC_HPCt2hg_1IQi,'all');
    PAC_HPCt2hg_means(2,ses)=mean(PAC_HPCt2hg_2IQi,'all');
    PAC_HPCt2hg_means(3,ses)=mean(PAC_HPCt2hg_3IQi,'all');
    PAC_HPCt2hg_means(4,ses)=mean(PAC_HPCt2hg_4IQi,'all');
    
    PAC_HPCt2hg_max(1,ses)=max(PAC_HPCt2hg_1IQi,[],'all');
    PAC_HPCt2hg_max(2,ses)=max(PAC_HPCt2hg_2IQi,[],'all');
    PAC_HPCt2hg_max(3,ses)=max(PAC_HPCt2hg_3IQi,[],'all');
    PAC_HPCt2hg_max(4,ses)=max(PAC_HPCt2hg_4IQi,[],'all');
    
    Time_IQi_THPWR_all{1,ses}=Time_1IQi_THPWR;
    Time_IQi_THPWR_all{2,ses}=Time_2IQi_THPWR;
    Time_IQi_THPWR_all{3,ses}=Time_3IQi_THPWR;
    Time_IQi_THPWR_all{4,ses}=Time_4IQi_THPWR;

end

%% I check the average speed of the animal for all the 500ms epochs that theta power is in the 1st, 2nd, 3rd or 4th interquartile (IQ) interval. FOR EVERY SESSION
% This is a control analysis to examine whether the differential modulation of gamma amplitude from theta phase depending on theta power is second to the speed of the animal
% Since all the segments are 500ms just the start is enough to index from the speed matrix Velocity_cmpersec_500
Vel_500_allThIQis=cell(size(Time_IQi_THPWR_all)); %that has the velocities of the individual epochs for all theta power IQis. No need to export that now
Vel_500_allThIQis_mean=NaN(size(Time_IQi_THPWR_all)); %this has the mean velocity of all theta power IQis for all sessions. %THIS NEED TO BE EXPORTED
for ses=1:size(Time_IQi_THPWR_all,2)
    for iqi=1:size(Time_IQi_THPWR_all,1)
        for epoch=1:size(Time_IQi_THPWR_all{iqi,ses},1)
            Vel_500_allThIQis{iqi,ses}(epoch,1)=Velocity_cmpersec_500(ses,(Time_500ms_all(1,:)==Time_IQi_THPWR_all{iqi,ses}(epoch,1)));
        end
        Vel_500_allThIQis_mean(iqi,ses)=mean(Vel_500_allThIQis{iqi,ses});
    end
end



%% Checking recurrence in PAC on a session by session basis for all Theta power IQis

PAC_TH_LG_cor_sesbyses= cell(1,size(PAC_HPCt2lg_all,1)); %the contain the autocor matrices for all IQis
PAC_TH_HG_cor_sesbyses= cell(1,size(PAC_HPCt2hg_all,1));

for IQi=1:size(PAC_HPCt2lg_all,1)

    for corX=1:sessions
        for corY=1:sessions
        PAC_TH_LG_cor_sesbyses{IQi}(corY,corX)= corr2(PAC_HPCt2lg_all{IQi,corX},PAC_HPCt2lg_all{IQi,corY});
        PAC_TH_HG_cor_sesbyses{IQi}(corY,corX)= corr2(PAC_HPCt2hg_all{IQi,corX},PAC_HPCt2hg_all{IQi,corY});
        end   
    end
    
    PAC_TH_LG_cor_sesbyses{IQi}= flipud(PAC_TH_LG_cor_sesbyses{IQi});
    PAC_TH_HG_cor_sesbyses{IQi}= flipud(PAC_TH_HG_cor_sesbyses{IQi});
end


% This was just a test to see if I can narrow down the PAC matrix area I
% take into account when I check for reccurence over the duration of the
% experiment.
% PAC_TH_LG_cor_Hotspot= cell(1,size(PAC_HPCt2lg_all,1));
% PAC_TH_HG_cor_Hotspot= cell(1,size(PAC_HPCt2hg_all,1));
% for IQi=1:size(PAC_HPCt2lg_all,1)
% 
%     [PACHotspots]=find(PAC_HPCt2lg_all{IQi,1}>=(median(PAC_HPCt2lg_all{IQi,1}, 'all'))); %find the frequency areas with highest MI values during the fist session 
% 
%     for corX=1:sessions
%         for corY=1:sessions
%         PAC_TH_LG_cor_Hotspot{IQi}(corY,corX)= corr2(PAC_HPCt2lg_all{IQi,corX}(PACHotspots),PAC_HPCt2lg_all{IQi,corY}(PACHotspots));
%         PAC_TH_HG_cor_Hotspot{IQi}(corY,corX)= corr2(PAC_HPCt2hg_all{IQi,corX}(PACHotspots),PAC_HPCt2hg_all{IQi,corY}(PACHotspots));
%         end   
%     end
%     
%     PAC_TH_LG_cor_Hotspot{IQi}= flipud(PAC_TH_LG_cor_Hotspot{IQi});
%     PAC_TH_HG_cor_Hotspot{IQi}= flipud(PAC_TH_HG_cor_Hotspot{IQi});
% end


%%  Not used right now %%%%%%% Calculate PAC using two methods for every minute of the experiment but not stratify based on Theta power.

% %Calculating Phase-Amplitude coherence. The method used is the Modulation
% %without thresholding for theta power.
% %Index (Canolty et al 2006 Science)
% 
% PAC_HPCt2HPClg_all_MI=cell(1,sessions*SesDurMin);
% PAC_HPCt2HPChg_all_MI=cell(1,sessions*SesDurMin);
% PAC_HPCt2HPClg__means_MI=zeros(1,sessions*SesDurMin);
% PAC_HPCt2HPChg__means_MI=zeros(1,sessions*SesDurMin);
% n=1; %a counter for index
% for ses=1:sessions
%     %Here I'm looking at minbymin epochs 
%     for fsec=1:size(EGF_500ms_all{ses},2)/10:size(EGF_500ms_all{ses},2)
%     %this is the main function for PAC from (Onslow,2010). The "repeated
%     %trials" for each 1min bin are 120 500ms epochs. 
%     %PAC between theta and low gamma
%     [PAC_HPCt2HPClg,~, ~]=find_pac_shf(EGF_500ms_all{ses}(:,fsec:fsec+size(EGF_500ms_all{ses},2)/10-1),Fs,'mi',EGF_500ms_all{ses}(:,fsec:fsec+size(EGF_500ms_all{ses},2)/10-1),Lims{1,2}:0.5:Lims{2,1},Lims{1,4}:0.5:Lims{2,3},'n',0,7,200,0,0.05,'HPC-HPC','HPC','HPC');
%     %PAC between theta and high gamma
%     [PAC_HPCt2HPChg,~, ~]=find_pac_shf(EGF_500ms_all{ses}(:,fsec:fsec+size(EGF_500ms_all{ses},2)/10-1),Fs,'mi',EGF_500ms_all{ses}(:,fsec:fsec+size(EGF_500ms_all{ses},2)/10-1),Lims{1,2}:0.5:Lims{2,1},Lims{1,6}:0.5:Lims{2,5},'n',0,7,200,0,0.05,'HPC-HPC','HPC','HPC');
%     
%     PAC_HPCt2HPClg_all_MI{1,n}=PAC_HPCt2HPClg;
%     PAC_HPCt2HPClg__means_MI(1,n)=mean(PAC_HPCt2HPClg,'all');
%     PAC_HPCt2HPChg_all_MI{1,n}=PAC_HPCt2HPChg;
%     PAC_HPCt2HPChg__means_MI(1,n)=mean(PAC_HPCt2HPChg,'all');
%     n=n+1;
%     end
% end
% 
% %this is to explore the recurrence of PAC patterns across the recording sessions. It might be interesting but probably not very sensitive measure.    
% PAC_TH_LG_cor_minbymin_MI=zeros(sessions*SesDurMin,sessions*SesDurMin);
% PAC_TH_HG_cor_minbymin_MI= zeros(sessions*SesDurMin,sessions*SesDurMin);
% for corX=1:sessions*SesDurMin
%     for corY=1:sessions*SesDurMin
%         PAC_TH_LG_cor_minbymin_MI(corY,corX)= corr2(PAC_HPCt2HPClg_all_MI{corX},PAC_HPCt2HPClg_all_MI{corY});
%         PAC_TH_HG_cor_minbymin_MI(corY,corX)= corr2(PAC_HPCt2HPChg_all_MI{corX},PAC_HPCt2HPChg_all_MI{corY});
%     end
% end
% PAC_TH_LG_cor_minbymin_MI=flipud(PAC_TH_LG_cor_minbymin_MI);
% PAC_TH_HG_cor_minbymin_MI=flipud(PAC_TH_HG_cor_minbymin_MI);
%         

%--------------------------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Calculating Phase-Amplitude coherenc. The method used is the Envelope-to-Signal Correlation  (Bruns et al 2004 Int J Psychophysiol)
% 
% PAC_HPCt2HPClg_all_ESC=cell(1,sessions*SesDurMin);
% PAC_HPCt2HPChg_all_ESC=cell(1,sessions*SesDurMin);
% PAC_HPCt2HPClg__means_ESC=zeros(1,sessions*SesDurMin);
% PAC_HPCt2HPChg__means_ESC=zeros(1,sessions*SesDurMin);
% n=1; %a counter for index
% for ses=1:sessions
%     %Here I'm looking at minbymin epochs 
%     for fsec=1:size(EGF_500ms_all{ses},2)/10:size(EGF_500ms_all{ses},2)
%     %this is the main function for PAC from (Onslow,2010). The "repeated
%     %trials" for each 1min bin are 120 500ms epochs. 
%     %PAC between theta and low gamma
%     [PAC_HPCt2HPClg,~, ~]=find_pac_shf(EGF_500ms_all{ses}(:,fsec:fsec+size(EGF_500ms_all{ses},2)/10-1),Fs,'esc',EGF_500ms_all{ses}(:,fsec:fsec+size(EGF_500ms_all{ses},2)/10-1),Lims{1,2}:0.5:Lims{2,1},Lims{1,4}:0.5:Lims{2,3},'n',0,7,200,0,0.05,'HPC-HPC','HPC','HPC');
%     %PAC between theta and high gamma
%     [PAC_HPCt2HPChg,~, ~]=find_pac_shf(EGF_500ms_all{ses}(:,fsec:fsec+size(EGF_500ms_all{ses},2)/10-1),Fs,'esc',EGF_500ms_all{ses}(:,fsec:fsec+size(EGF_500ms_all{ses},2)/10-1),Lims{1,2}:0.5:Lims{2,1},Lims{1,6}:0.5:Lims{2,5},'n',0,7,200,0,0.05,'HPC-HPC','HPC','HPC');
%     
%     PAC_HPCt2HPClg_all_ESC{1,n}=PAC_HPCt2HPClg;
%     PAC_HPCt2HPClg__means_ESC(1,n)=mean(PAC_HPCt2HPClg,'all');
%     PAC_HPCt2HPChg_all_ESC{1,n}=PAC_HPCt2HPChg;
%     PAC_HPCt2HPChg__means_ESC(1,n)=mean(PAC_HPCt2HPChg,'all');
%     n=n+1;
%     end
% end
% 
% %this is to explore the recurrence of PAC patterns across the recording sessions. It might be interesting but probably not very sensitive measure.    
% PAC_TH_LG_cor_minbymin_ESC=zeros(sessions*SesDurMin,sessions*SesDurMin);
% PAC_TH_HG_cor_minbymin_ESC= zeros(sessions*SesDurMin,sessions*SesDurMin);
% for corX=1:sessions*SesDurMin
%     for corY=1:sessions*SesDurMin
%         PAC_TH_LG_cor_minbymin_ESC(corY,corX)= corr2(PAC_HPCt2HPClg_all_ESC{corX},PAC_HPCt2HPClg_all_ESC{corY});
%         PAC_TH_HG_cor_minbymin_ESC(corY,corX)= corr2(PAC_HPCt2HPChg_all_ESC{corX},PAC_HPCt2HPChg_all_ESC{corY});
%     end
% end
% PAC_TH_LG_cor_minbymin_ESC=flipud(PAC_TH_LG_cor_minbymin_ESC);
% PAC_TH_HG_cor_minbymin_ESC=flipud(PAC_TH_HG_cor_minbymin_ESC);

%% -Here I calculate instateneous Amp and Phase that I'll be using in the spike phase locking section------------------------------------------------------------------------------%%%%%%%%%%%%%%%%%%%


% Band pass filter. I use the frequency limits set at the start of this
% script and apply them to the raw LFP for as many times as the frequency
% bands I'm interested in. In this case I'm interested in theta low-gamma
% and high gamma.
%this is a cell array with 3 rows 6 columns. Each row correspond to an LFP band and each column corresponds to a session. Top row is theta, second row is slow gamma and third row is high gamma.

EGF_bandpass_all= cell(length(Lims)/2, sessions); 
for ses=1:sessions
    indx=1;
    for i=1:2:length(Lims)-1
    EGF_bandpass = fftbandpass(EGF_all_uV{ses}(1:Fs*60*SesDurMin),Fs, Lims{1,i}, Lims{1,i+1}, Lims{2,i}, Lims{2,i+1}); % Theta, low-gamma and high-gamma bands
    EGF_bandpass_all{indx,ses}=EGF_bandpass'; %I invert just because I like column vectors
    indx=indx+1;
    end
end

% Calculate eeg instantaneous phase based on TROUGHs with Hilbert transform
disp('Computing instantaneous Phase & Amplitute with Hilbert transform for all sessions and oscillation bands')
Phase_all=cell(length(Lims)/2,sessions);
Amp_all=cell(length(Lims)/2, sessions);
for ses=1:sessions
    for i=1:length(Lims)/2
        [phase,Amp] = Osciphase(EGF_bandpass_all{i,ses});
        Phase_all{i,ses}=phase;
        Amp_all{i,ses}=Amp;
    end
end

% Instantaneous power for the oscillation bands of interest. This is a
% temporary solution; wavelet analysis could be better.
Instpower_all=cell(length(Lims)/2,sessions);
for ses=1:sessions
    for i=1:length(Lims)/2
        instpower = Amp_all{i, ses}.^ 2;
        Instpower_all{i,ses}=instpower;
    end
end

disp('Computing Theta and Gamma Phase assymetries in all sessions')
% This is related to Trimper et 2014 Hippocampus (Fig6) and Kitanashi et al
% 2015 Neuron (Fig2H. Here I use the KItanishi et al 30deg bins (uniform
% distribution 8.3%)
Th_Ph_bins=0:30:360;
Phase_portion_all=cell(3,6); Phase_portion_all(:)={NaN(1,length(Th_Ph_bins)-1)};
for ses=1:sessions
    for band=1:length(Lims)/2    
        for ph=1:length(Th_Ph_bins)-1
        %this calculates the porportion of a certain phase of an
        %oscillation compared to the entirety of the signal. For 90deg bin
        %that should be .25 (Trimper et al 2014) is the wave was
        %symmetrical
        Phase_portion_all{band,ses}(ph)= length(Phase_all{band,ses}(Phase_all{band,ses}>=Th_Ph_bins(ph) & Phase_all{band,ses}<Th_Ph_bins(ph+1)))/length(Phase_all{band,ses}); 
        end
    end
end

      
%These cell arrays have the same structure as and every cells has the same vector dimensions. 
%I'll use that to find the mean amp of gamma (low and high) for every phase bin in the different power ranges of theta
disp('Computing Theta Power interquartile (IQ) intervals for all sessions')
Theta_Ph_IQi_all=cell(4,sessions); %It contains the phase of theta for different theta powers. Each row is a different interquartile (IQ) interval of theta power. Each column is a session.  
LG_Amp_IQi_all=cell(4, sessions); %It contains the amplitudes of low gamma for different theta power distributions. 
HG_Amp_IQi_all=cell(4, sessions);
for ses=1:sessions

    [q1, q2, q3] = findiquartile(Instpower_all{1,ses}); %finds the thresholds for 1st, 2nd, 3rd and 4th interquartile (IQ) interval. This is based on Tort et al 2009 PNAS (Fig 4A more specifically)

    Phase_1IQi=Phase_all{1,ses}(Instpower_all{1,ses}<q1);
    LG_Amp_1IQi=Amp_all{2,ses}(Instpower_all{1,ses}<q1);
    HG_Amp_1IQi=Amp_all{3,ses}(Instpower_all{1,ses}<q1);
    Phase_2IQi=Phase_all{1,ses}(Instpower_all{1,ses}>=q1 & Instpower_all{1,ses}<q2);
    LG_Amp_2IQi=Amp_all{2,ses}(Instpower_all{1,ses}>=q1 & Instpower_all{1,ses}<q2);
    HG_Amp_2IQi=Amp_all{3,ses}(Instpower_all{1,ses}>=q1 & Instpower_all{1,ses}<q2);
    Phase_3IQi=Phase_all{1,ses}(Instpower_all{1,ses}>q2 & Instpower_all{1,ses}<q3);
    LG_Amp_3IQi=Amp_all{2,ses}(Instpower_all{1,ses}>q2 & Instpower_all{1,ses}<q3);
    HG_Amp_3IQi=Amp_all{3,ses}(Instpower_all{1,ses}>q2 & Instpower_all{1,ses}<q3);
    Phase_4IQi=Phase_all{1,ses}(Instpower_all{1,ses}>=q3);
    LG_Amp_4IQi=Amp_all{2,ses}(Instpower_all{1,ses}>=q3);
    HG_Amp_4IQi=Amp_all{3,ses}(Instpower_all{1,ses}>=q3);
    %catching all the outputs
    Theta_Ph_IQi_all{1,ses}=Phase_1IQi; Theta_Ph_IQi_all{2,ses}=Phase_2IQi; Theta_Ph_IQi_all{3,ses}=Phase_3IQi; Theta_Ph_IQi_all{4,ses}=Phase_4IQi;
    LG_Amp_IQi_all{1,ses}=LG_Amp_1IQi; LG_Amp_IQi_all{2,ses}=LG_Amp_2IQi; LG_Amp_IQi_all{3,ses}=LG_Amp_3IQi; LG_Amp_IQi_all{4,ses}=LG_Amp_4IQi;
    HG_Amp_IQi_all{1,ses}=HG_Amp_1IQi; HG_Amp_IQi_all{2,ses}=HG_Amp_2IQi; HG_Amp_IQi_all{3,ses}=HG_Amp_3IQi; HG_Amp_IQi_all{4,ses}=HG_Amp_4IQi;
    
end


disp('Calculating mean gamma amplitude for theta phases in all theta power interquartile (IQ) intervals')
Phase_bins=0:20:360;
Mean_GAmp_perTHphasebin_all=cell(4, sessions); %4 interquartile (IQ) intervals and 6 sessions
Mean_GAmp_perTHphasebin_all(:)={NaN(2,length(Phase_bins)-1)}; %2 bands of gamma and 18 theta phase bins
for ses=1:sessions %for all sessions
    for IQi=1:size(Theta_Ph_IQi_all,1) %for all theta power interquartile (IQ) interval
        for ph=1:length(Phase_bins)-1 %for all the different theta phase bins
            
            mean_LGAmp_perTHphasebin=mean(LG_Amp_IQi_all{IQi,ses}(Theta_Ph_IQi_all{IQi,ses}>=Phase_bins(ph) & Theta_Ph_IQi_all{IQi,ses}<Phase_bins(ph+1)));
            mean_HGAmp_perTHphasebin=mean(HG_Amp_IQi_all{IQi,ses}(Theta_Ph_IQi_all{IQi,ses}>=Phase_bins(ph) & Theta_Ph_IQi_all{IQi,ses}<Phase_bins(ph+1)));
            
            Mean_GAmp_perTHphasebin_all{IQi,ses}(1,ph)=mean_LGAmp_perTHphasebin; %catch the low gamma amplitude
            Mean_GAmp_perTHphasebin_all{IQi,ses}(1,ph+length(Phase_bins)-1)=mean_LGAmp_perTHphasebin; %the only reason I'm doing that is to get an overall 720deg plot and see the modulation better
            Mean_GAmp_perTHphasebin_all{IQi,ses}(2,ph)=mean_HGAmp_perTHphasebin; %catch the high gamma amplitude
            Mean_GAmp_perTHphasebin_all{IQi,ses}(2,ph+length(Phase_bins)-1)=mean_HGAmp_perTHphasebin;
        end
    end
end



%Maybe I should test for circular uniformity_V-test (expected angle) 



% Detect peak timing & temporal length of gamma events [sec]
%the name of the function is misleading since I'm detecting events for all
%oscillatory bands of interest. Not sure where will I use the theta events though.
%The slow and fast gamma events will be used in the main function on the final stage where spikes and LFP will come together.
%Kitanishi et al 2015 Neuron only use spike events during strong gamma to
%calculate locking
disp('Detecting theta and gamma events for all sessions')
Peaks_all=cell(length(Lims)/2,sessions);
Length_all =zeros(length(Lims)/2,sessions);
for ses=1:sessions
    for i=1:length(Lims)/2
        [peaks, tlength] = gammaEvent(Instpower_all{i,ses},Fs,wLength);
        Peaks_all{i,ses}=peaks;
        Length_all(i,ses)=tlength;
    end
end


%-------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section Three. Importing spike timings for the clusters indicated in the Clu2use xl file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Reading cluster identities from .clu files. These are all the clusters that were detected in this session

disp('Reading cluster identities from .clu files...');
CluID_all = cell(1,8); %a 1x8 cell array with the cluster indentities for all 8 tetrodes
for i = 1:length(elecs)
        current_e = elecs(i);
clufile = ['merge.clu.' num2str(current_e)];
fid_clu = fopen(clufile,'r');
disp(['Getting cluster identities for tetrode ', num2str(i)]);
nclus = sscanf(fgetl(fid_clu),'%d');
cluID_elec = fscanf(fid_clu,'%f',[1 Inf]);
CluID_all{i}=cluID_elec;
fclose(fid_clu); %close the file at the end of each cycle to make the process faster
end

%Getting all the spike features from this session.

disp('Reading spike features from .fet files...this might take a while');
Feat_all = cell(1,8); %a 1x8 cell array with the spike features for all 8 tetrodes
for k = 1:length(elecs)
        current_e = elecs(k);
fetfile = ['merge.fet.' num2str(current_e)];
fid_fet = fopen(fetfile, 'r');
disp(['Getting features for tetrode ', num2str(k)]);
nfeatures = sscanf(fgetl(fid_fet),'%d');
feat_elec = fscanf(fid_fet,'%f',[nfeatures Inf]);
Feat_all{k}=feat_elec;
fclose(fid_fet);
end

% Gettinh the usable clusters from the manually generated xl file. These
% clusters were also used for the ensemble coordination analysis

disp('Choosing usable clusters from xl file...');
Clu2useID = xlsread('Clu2use.xlsx'); %open the file which contains the usable clusters-needs to be created manually still
Tetr= unique(Clu2useID(:,1)); %find the tetrodes that had good clusters
disp([num2str(length(Clu2useID)), ' good clusters in total in this session']);
GoodCluID_all = cell(1,8); %a 1x8 cell array with the good cluster indentities for all 8 tetrodes
for n=1:length(Tetr) %finding the usable clusters for each tetrode and creates a list for each tetrode
I=find(Clu2useID(:,1)==Tetr(n));
goodclusID = Clu2useID(I(1):I(end),2);
GoodCluID_all{Tetr(n)}= goodclusID;
end

%Using the information above to index the cells arrays with the info for
%all clusters and get the index of them

Tetr_goodclus=find(~cellfun(@isempty,GoodCluID_all)); %which tetrodes have good clusters
GoodCluindx_all = cell(1,length(Clu2useID)); %a 1xn cell array with all the good cluster indentities from all tetrodes
h=0; % I need a counter that will count the total cluster number so I can index the cell array correctly
for f=1:length(Tetr_goodclus) %going through the tetrodes with the good clusters
    for c=1:length(GoodCluID_all{Tetr_goodclus(f)})
        gcInx=find(CluID_all{Tetr_goodclus(f)}==GoodCluID_all{Tetr_goodclus(f)}(c)); %finding the index of the good clusters in the vectors from .clu files
        h=h+1;
        GoodCluindx_all{h}=gcInx;
    end
end

% Use the indexes found above and get the spike timings that correspond to
% the "good" clusters from the 24 row of the corresponding feature matrix

TtrID = Clu2useID(:,1);
T4GoodClu=cell(1,length(TtrID)); %a cell array with the same dimensions as GoodCluindx_all
for ii=1:length(TtrID)
    GT=Feat_all{TtrID(ii)}(24,:);
    t4goodclu=GT(GoodCluindx_all{ii});
    T4GoodClu{ii}=t4goodclu;
end


%Remember!! The spike timings are for the combination of all sessions while the LFP
% data times are separate.


%I'm using the T_segs now to correct the spike timings and separate them in
%sessions. Keep in mind that there is a padding of 10sec for every session.
%While the duration of individual sessions is 600.96sec the difference
%between the starts of each session is 610.96sec.
T4GoodClu_sep =cell(length(T4GoodClu), sessions); %this is a cell array that contains the spk times of each cluster (row) for all sessions(columns).
NoSpk_cluses=[]; %This is to catch the number cluster that did not spike in a session. That I don't expect this to happen often therefore I am not pre allocating space_Also I don't know.
for clu=1:length(T4GoodClu)
    for tp=1:length(T_segs)-1
        inds = not(abs(sign(sign(T_segs(tp) - T4GoodClu{clu}) + sign(T_segs(tp+1) - T4GoodClu{clu})))); %finding the index of of spike times in the T4GoodClu that are between the start and the end of the segment
        if sum(inds)==0
            idx=1; %start a counter for every time the above condition is met so I can index a matrix and keep the session number(s) a cluster did not spike in.
            disp(['Cluster ', num2str(clu), ' did not have any spikes in session ',num2str(tp)]); %Sometimes a cluster hasn't fired in all sessions.
            NoSpk_cluses(idx,1)=clu; %This might be used later so the appropriate cell will be skipped
            NoSpk_cluses(idx,2)=tp;
            idx=idx+1;
        end
        spk_ses=T4GoodClu{clu}(inds); %getting the spike timings in that session
        spk_ses=spk_ses-T_segs(tp); %correcting for the start of each session which includes the padding
        T4GoodClu_sep{clu,tp}=spk_ses;
    end
end



%The output of this loop 'T4GoodClu_sep' contains the timestamps for all the spikes of the good quality clusters with corrected time.
%I'll check the speed of the animal and exclude spikes that happened during
%periods of imobility. I'll then see what is the instateneous phase for each
%oscillation fo interest (theta, slow gamma, fast gamma) when the spikes happened and check if there is phase locking.

%% Section Four %%%% Bringing spikes and LFP together

%I now take the spike timestamps from T4GoodClu_sep and find what is the phase of theta, slow gamma and fast gamma at the time each
%spike happened based on the corresponding Phase_all cell. The output is three cell arrays each for every oscillatory
%band. They should have the same dimensions as T4GoodClu_sep and
%all their cells have the same elements with the equivalent cells in
%T4GoodClu_sep.


%1/Fs is the time increment based on the sampling frequency
% Because of the frequenccy range the first data point will be 1/Fs sec after
% the start of the recording. In total we have 601/(1/Fs) sampling points
% that 2884800 exactly the length of EGF
%Creation of a timevector that will help me index the spike timestamps.
%All elements in Phase_all should have the same length therefore I don;t
%have to create it 18 times
disp('Finding oscillation phases for every spike....this might take a while');

Tvec = 1/Fs + (0:length(Phase_all{1,1})-1)*(1/Fs);
%this appoach with multiple for loops is very very slow...I need to revisit!
SpkPhaseT_all=cell(size(T4GoodClu_sep));
SpkPhaseS_all=cell(size(T4GoodClu_sep));
SpkPhaseF_all=cell(size(T4GoodClu_sep));
for ses=1:sessions %for all sessions
    for clu=1:size(T4GoodClu_sep, 1) %for all clusters
        if isempty(T4GoodClu_sep{clu,ses}) %this is to skip to the next iteration of the loop if a cluster doesn't have spike during a specific session
            continue
        end
        disp(['Finding oscilation phases for cluster ', num2str(clu), ' in session ',num2str(ses)]);
        for spkT=1:length(T4GoodClu_sep{clu,ses})
            [~, ix] = min(abs(Tvec-(T4GoodClu_sep{clu, ses}(spkT))/100000)); %find the LFP sampling time that corresponds (or it's the closest) to this spike timestamp
            %find the phase for all different bands of oscillation
            SpkPhaseT= Phase_all{1,ses}(ix); %for theta
            SpkPhaseT_all{clu,ses}(spkT)=SpkPhaseT; %catch it in the right cell array
            SpkPhaseS= Phase_all{2,ses}(ix); %for slow gamma
            SpkPhaseS_all{clu,ses}(spkT)=SpkPhaseS;
            SpkPhaseF= Phase_all{3,ses}(ix); %for fast gamma
            SpkPhaseF_all{clu,ses}(spkT)=SpkPhaseF;
        end
    end
end


%% Now I'll get the spikes that happend during gamma events. For theta is the entire session
%With that I can calculate the 

%from Laura Colgin 2009 Nature
% Cells were considered to be phase
% -locked to slow or fast gamma if their phase distribution 
% differed significantly from uniform (P < 0.0
% 5, Rayleigh test).

%I do screen for gamma events (strong gamma) though so I don't think I should screen for speed as well as gamma event 
SPKPhaseT=cell(size(T4GoodClu_sep));
SPKPhaseS=cell(size(T4GoodClu_sep));
SPKPhaseF=cell(size(T4GoodClu_sep));
Rate_All=cell(size(T4GoodClu_sep)); %this has the firing rate of each cluster in each session during the entire session (1st row), slow gamma events only (2nd row), and fast gamma event only (3rd row)
SPKRadPhaseT=cell(size(T4GoodClu_sep)); %these should have identical identical dimensions with SPKPhase cell arrays but values in radians instead of degrees
SPKRadPhaseS=cell(size(T4GoodClu_sep));
SPKRadPhaseF=cell(size(T4GoodClu_sep));
ResulV_all=cell(size(T4GoodClu_sep)); %Resultant vector lengths for every cluster and every session, for theta, slow gamma(2nd row), and fast gamma(3rd row)
MeanAng_all=cell(size(T4GoodClu_sep)); %Same dimensions and strusture as ResulV_all. It contains the mean resultant phase
Raylp_all=cell(size(T4GoodClu_sep)); %p values for Rayleigh test
NormSpk_hist_all=cell(size(T4GoodClu_sep)); %this has normalised spkcounts for different phases (30deg bins) of theta(first row of each cell), slow gamma(second row of each cell) and fast gamma(third row of each cell). For every cluster(rows) and every session(columns). Sessions for clusters that spiked less than 36 spikes will be excluded and be empty
for clu=1:size(T4GoodClu_sep, 1)%for all clusters
    for ses=1:sessions %for all sessions
        if length(T4GoodClu_sep{clu,ses})<36 %this is to skip to the next iteration of the loop if a cluster has less than 36 spikes (3 times the number of phase bins) during a specific session
            continue
        end
        disp(['Finding spikes for cluster ', num2str(clu), ' that happened during gamma evens in session ',num2str(ses)]);
        idxS = selectSpikes(T4GoodClu_sep{clu,ses}./100000, Peaks_all{2,ses}, wLength);
        idxF = selectSpikes(T4GoodClu_sep{clu,ses}./100000, Peaks_all{3,ses}, wLength);
        SPKPhaseT{clu,ses} = SpkPhaseT_all{clu,ses}; %No spike exclusion for Theta. But session for clusters that spiked less than 36 spikes will be excluded.
        SPKPhaseS{clu,ses} = SpkPhaseS_all{clu,ses}(logical(idxS)); %the slow gamma phases of spikes that happened during slow gamma events 
        SPKPhaseF{clu,ses} = SpkPhaseF_all{clu,ses}(logical(idxF)); %the high gamma phases of spikes that happened during high gamma events 
        
        % "Firing rate" of each cluster that is calculated based on the
        % spikes that happened during the events divided by the total
        % length of these events. For theta there's no excusion so the
        % duration is the duration of the trial
        Rate_All{clu,ses}(1,:) = length(SPKPhaseT{clu,ses})/SesDurMin*60; %that is the duration of the trial 600sec
        Rate_All{clu,ses}(2,:) = length(SPKPhaseS{clu,ses})/Length_all(2,ses); %these I get from the gamma event
        Rate_All{clu,ses}(3,:) = length(SPKPhaseF{clu,ses})/Length_all(3,ses);
        
        % Degree --> Radian spike phase. This is from the circular statistics toolbox. Nothing complicated though.
        SPKRadPhaseT{clu,ses} = circ_ang2rad(SPKPhaseT{clu,ses});
        SPKRadPhaseS{clu,ses} = circ_ang2rad(SPKPhaseS{clu,ses});
        SPKRadPhaseF{clu,ses} = circ_ang2rad(SPKPhaseF{clu,ses});
        
        % Resultant vector length. That's again from the circular
        % statistics toolbox. I will use these values to plot as a bar
        % graph the average of all the cell from all the animals 
        ResulV_all{clu,ses}(1,:) = circ_r(SPKRadPhaseT{clu,ses},[], [], 2);
        ResulV_all{clu,ses}(2,:) = circ_r(SPKRadPhaseS{clu,ses},[], [], 2);
        ResulV_all{clu,ses}(3,:) = circ_r(SPKRadPhaseF{clu,ses},[], [], 2);
        
        % Resultant vector phase [deg]. Notice that I have to calculate the
        % mean in radians and then covert back to degrees (more intuitive
        % for us stupid biologists)
        
        MeanAng_all{clu,ses}(1,:) = circ_rad2ang(circ_mean(SPKRadPhaseT{clu,ses}, [], 2));
        MeanAng_all{clu,ses}(2,:) = circ_rad2ang(circ_mean(SPKRadPhaseS{clu,ses}, [], 2));
        MeanAng_all{clu,ses}(3,:) = circ_rad2ang(circ_mean(SPKRadPhaseF{clu,ses}, [], 2));

        % Rayleigh test (circ_rtest). That's for each unit and each
        % session. Sessions in which clusters fired less than 36 spikes are
        % excluded
        
        Raylp_all{clu,ses}(1,:) = circ_rtest(SPKRadPhaseT{clu,ses});
        Raylp_all{clu,ses}(2,:) = circ_rtest(SPKRadPhaseS{clu,ses});
        Raylp_all{clu,ses}(3,:) = circ_rtest(SPKRadPhaseF{clu,ses});
             
        % Organise spike phases in a histogram and normalize spike count based on the number of spikes in the session/period of storng gamma. 
        x=15:30:345; %In hist 'x' is a vector of evenly spaced values, then hist uses the values as the bin centers.

        
        %I need to change that to histogram OR find a different way to bin
        %the normalised spikes
        nT = hist(SPKPhaseT{clu,ses},x); %this 
        nS = hist(SPKPhaseS{clu,ses},x);
        nF = hist(SPKPhaseF{clu,ses},x);
        NormSpk_hist_all{clu,ses}(1,:) = nT./length(SPKPhaseT{clu,ses}); %normalised to the number of spikes that happened in the sessions(no theta power screening)
        NormSpk_hist_all{clu,ses}(2,:) = nS./length(SPKPhaseS{clu,ses}); %normalised to the number of spikes that happened during slow gamma events
        NormSpk_hist_all{clu,ses}(3,:) = nF./length(SPKPhaseF{clu,ses}); %normalised to the number of spikes that happened during fast gamma events
        
    end
end

%% Exporting outputs

%First for the power spectral density analysis
mkdir LFPAnalysisOutput %creates the folder in which the figures and xl files will be stored

baseFileName = 'PWRanalysis.xls'; %creating file that contains the data for the power septral density analysis. Section Two-line 317 onwards
fullFileName = fullfile('LFPAnalysisOutput', baseFileName);
csvwrite(fullfile('LFPAnalysisOutput','Vel_500ms.csv'),Velocity_cmpersec_500);
xlswrite(fullFileName ,horzcat(f_spect',spectrogram_output.spectro.mean),'MeanSpectro');
xlswrite(fullFileName ,horzcat(f_spect',spectrogram_output.db_conv.mean),'MeanSdB_conv');
xlswrite(fullFileName ,horzcat(f_spect',spectrogram_output.spectro_norm.mean),'MeanSpectro_norm');


%this is for individual bands theta, slow gamma and fast gamma. On a minute
%by minute 
for ses=1:sessions
    xlswrite(fullFileName ,vertcat(frequency_bands{1,ses}.mean_theta(1:SesDurMin),frequency_bands{1,ses}.mean_slow_gamma(1:SesDurMin),frequency_bands{1,ses}.mean_fast_gamma(1:SesDurMin)),['PWR_fbands_Session', num2str(ses)]);
end

%this create a sheet with the normalised power for the three bands of interest and the LG/HG ratio in
%different velocity bins.  Data from line434
v=0:40;
vm=0.5:39.5; %that has the middle of the velocity bins so I can plot
n=1;
for ses=1:sessions
    meanpwr=NaN(4,length(v)-1); %create the matrix that will be exported
    for band=1:1:length(Lims)/2
        
        for vbin=1:length(v)-1
            
            meanpwr(band,vbin)=mean(PWR_bands_500_nooutliers{band,ses}(Vel_4_PWR_bands{band,ses}>=v(vbin) & Vel_4_PWR_bands{band,ses}<v(vbin+1)));
            meanpwr(4,vbin)=mean(PWR_LG2HG_500_ratio_nonoutliers{ses}(Vel_4_LG2HG_500_ratio{ses}>=v(vbin) & Vel_4_LG2HG_500_ratio{ses}<v(vbin+1)));
      
        end
  
    end
    
    xlswrite(fullFileName ,vertcat(vm,meanpwr),'PWRsvsVel',['A', num2str(n)]);
    n=n+6; 
end

%Now for the PAC (theta-gamma) related analysis.
baseFileName = 'PACrelatedanalysis.xls'; %creating file that contains the data for the Phase Amplitude coherence related analyses. Section Two-line 450 onwards
fullFileName = fullfile('LFPAnalysisOutput', baseFileName);

%this have the PACmeans and max for all 4 theta IQis and all sessions. Also
%a control analysis to see if the speed is the modulator of potential
%effects. This is related to Tort et al 2009 PNAS Figure 4B
xlswrite(fullFileName ,PAC_HPCt2lg_means,'PACmeans','A1');
xlswrite(fullFileName ,PAC_HPCt2hg_means,'PACmeans','A6');
xlswrite(fullFileName ,PAC_HPCt2lg_max,'PACmax','A1');
xlswrite(fullFileName ,PAC_HPCt2hg_max,'PACmax','A6');
xlswrite(fullFileName ,Vel_500_allThIQis_mean,'Vel_mean4ThetaIQis');


Ph_bin_middle=15:30:345;
%This is for the oscillation phase asymmetries_related to Kitanishi et al
%Figure 2H
xlswrite(fullFileName ,Ph_bin_middle,'ThetaPhasesymmetry','A1');
xlswrite(fullFileName ,Ph_bin_middle,'LgammaPhasesymmetry','A1');
xlswrite(fullFileName ,Ph_bin_middle,'HgammaPhasesymmetry','A1');
for ses=1:sessions
    
    xlswrite(fullFileName ,Phase_portion_all{1,ses},'ThetaPhasesymmetry',['A', num2str(ses+1)]);
    xlswrite(fullFileName ,Phase_portion_all{2,ses},'LgammaPhasesymmetry',['A', num2str(ses+1)]);
    xlswrite(fullFileName ,Phase_portion_all{3,ses},'HgammaPhasesymmetry',['A', num2str(ses+1)]);
end

%This is again for the modulation of gamma amplitude by theta phase. It's
%related to Tort et al 2009 PNAS Figure 4A
Phase_bins_middle=10:20:710;

xlswrite(fullFileName,Phase_bins_middle,'MeanGAmpThphase','A1'); %Thats to indicate the bins
n=2;
%with this loop the structure of the sheet will be lowgamma, highgamma mean amplitude in 4 theta power IQis for session 1 and then again for session 2 until session 6 
for ses=1:size(Mean_GAmp_perTHphasebin_all,2) 
    for iqi=size(Mean_GAmp_perTHphasebin_all,1)
        
        xlswrite(fullFileName ,Mean_GAmp_perTHphasebin_all{iqi,ses},'MeanGAmpThphase',['A', num2str(n)]);
        n=n+3; 
        
    end
end

%With this I collect the info I need so I can plot things like Kitanishi et al 2015 Neuron Figure 3B,D,E,F          
for ses=1:size(ResulV_all,2)
    n=1;
    for clu=1:size(ResulV_all,1)
        if isempty(ResulV_all{clu,ses})
            continue
        end
        
        xlswrite(fullFileName ,ResulV_all{clu,ses},['SpkPhaseLock_Session_', num2str(ses)],['A', num2str(n)]);
        xlswrite(fullFileName ,MeanAng_all{clu,ses},['SpkPhaseLock_Session_', num2str(ses)],['B', num2str(n)]);
        xlswrite(fullFileName ,Raylp_all{clu,ses},['SpkPhaseLock_Session_', num2str(ses)],['C', num2str(n)]);
        
        n=n+4;
    end
end


%With this I collect the info I need so I can plot things like Kitanishi et
%al 2015 Neuron Figure 3A. 

for ses=1:size(NormSpk_hist_all,2)
    xlswrite(fullFileName ,[x,x],['NormSpkCount_Session_', num2str(ses)],'A1');
    n=2;
    for clu=1:size(NormSpk_hist_all,1)
        if isempty(NormSpk_hist_all{clu,ses})
            continue
        end
        
        xlswrite(fullFileName ,[NormSpk_hist_all{clu,ses},NormSpk_hist_all{clu,ses}],['NormSpkCount_Session_', num2str(ses)],['A', num2str(n)]);
        
        n=n+4;
    end
end

%here I save the variabls that I didn't export in case Emma and Peter want
%something else in the analysis:

save([experiment,'extravariables','.mat'],'spectrogram_output','frequency_bands', 'theta_spect','slow_gamma_spect' ,'fast_gamma_spect','PAC_HPCt2lg_all','PAC_HPCt2hg_all','PAC_TH_LG_cor_sesbyses','PAC_TH_HG_cor_sesbyses','Rate_All','-v7.3','-nocompression')  

%% Exporting figures

mkdir LFPAnalysisFigures

for ses=1:size(StrongLFP_bandpass,2)
   
    fig=figure;set(gcf,'Visible', 'off');       
    subplot(4,1,1);    
    plot(StrongLFP_bandpass{1, ses});     
    subplot(4,1,2);    
    plot(StrongLFP_bandpass{2, ses});
    subplot(4,1,3);    
    plot(StrongLFP_bandpass{3, ses});
    subplot(4,1,4);    
    plot(StrongLFP_bandpass{4, ses});
    
    saveas(gcf,fullfile('LFPAnalysisFigures',['Example traces for Session_',num2str(ses)]),'svg');
    close
end



disp('FINISHED!')
toc;
%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%These are the custom functions that will be used in this analysis script.
%A few will be from the toolboxes FMAT, chronux, PAC and fieldtrip. I'll
%comment to indicate usage of toolbox functions

%% startup function for figures

function startup
%Figure defaults
%axes
set(0, 'DefaultTextFontSize', 14);
set(0, 'DefaultAxesFontSize', 16);
set(0, 'DefaultUicontrolFontSize', 8);
set(0,'DefaultAxesTickDir', 'out')
set(0,'DefaultAxesLineWidth', 1.5)

%data
set(0,'DefaultLineLineWidth',1)
co = [0    0    1
      1    0    0
    0.3  0.3  0.3
    0.2    1  0.2
    0      1    1
    0      0    0];
set(0,'DefaultAxesColorOrder',co)
set(0,'DefaultAxesAmbientLightColor',[1 1 1])
set(0, 'DefaultAxesBox', 'off')
end


%% readEGF
% Reads high sampled eeg data and returns it in the array EEG together with the
% sampling rate Fs in Hz

function [EEG,Fs] = readEGF(datafile)

% Open the file for reading
fid = fopen(datafile,'r');
% Skip some lines of the header
for ii = 1:8
   textstring = fgetl(fid);
end

% Read out the sampling rate
Fs = sscanf(textstring(13:end-3),'%f');
% Read out the number of bytes per sample
textstring = fgetl(fid);
bytesPerSample = sscanf(textstring(18:end),'%f');
% Skip some more lines of the header
textstring = fgetl(fid);
% Read out the number of samples in the file
nosamples = sscanf(textstring(17:end),'%u');
% Go to the start of data (first byte after the data_start marker)
fseek(fid,10,0);

% Read data according to the number of bytes used per sample
if bytesPerSample == 1
    if isempty(nosamples) %for some sessions the information for the number of samples is not there or it's unobtainable
        EEG = fread(fid,'int8');
        disp('Number of samples information could not be found, all egf data points were loaded for this session');
    else    
        EEG = fread(fid,nosamples,'int8');
    end
else
    if isempty(nosamples)
        EEG = fread(fid,'int16');
        disp('Number of samples information could not be found, all egf data points were loaded for this session');
    else
        EEG = fread(fid,nosamples,'int16');
    end
end
% Close the file
fclose(fid);
end

%% EegMode Takuma Kitanishi 101025

function [eegmode, gain] = EegMode(eeg)

% get eeg info. This is taken from the set file. Files were recorded in two
% versions of Axona recording systems therefore I have a bit at the start
% to check what system the data are coming from.
% Normally eeg mode should be 1 (Ref).
% If eeg mode is 3 (-Ref), eeg sign should be inverted before later
% analysis.

daqversion=0; %old system, 1= new system
fid2=fopen(eeg,'r');
tekst=fscanf(fid2, '%c');

refindeks=strfind(tekst, 'ref_'); %this a general comment for strfind function. The syntax at some point chenged and the bigger string should be the first input not the second.
if length(refindeks)>1
     disp('Recording was made with new system');
    daqversion=1; % new system
end

indeks=strfind(tekst, 'saveEEG');
antalleeg=0;
for i=1:length(indeks)
%     disp(tekst(indeks(i):indeks(i)+13));
    if tekst(indeks(i)+13)=='1'
        antalleeg=antalleeg+1;
    elseif tekst(indeks(i)+13)==0
        disp(strcat,'EEG number ',num2str(i),' was not saved.');
    end
end
% disp(strcat(num2str(antalleeg),' EEG files was saved.'));

kindeks=strfind(tekst, 'EEG_ch');
if daqversion==1; kindeks=kindeks(1:end-2); end

j=1;
for i=1:antalleeg
   %    disp(tekst(kindeks(j):kindeks(j)+10));

   %
   kanal=tekst(kindeks(j)+9:kindeks(j)+10);
   eegch(i)=str2double(kanal);
   kanalindeks=str2double(kanal)-1;

   % Search eeg mode ('ref etc...')
   searchstring=strcat('mode_ch_',num2str(kanalindeks));
   modeind=strfind(tekst, searchstring);
   tallmode=str2double(tekst(modeind(1)+10:modeind(1)+12));
   eegmode(i)=tallmode;

   % Search gain
   searchstring2 = strcat('gain_ch_',num2str(kanalindeks));
   gainind = strfind(tekst, searchstring2);
   gainindst = gainind(1)+10;   % gain read start ind
   gainindend = gainindst+1; % gain read end ind

   while 1
       if tekst(gainindend) == sprintf('\r') || tekst(gainindend) == sprintf('\n')
           gainindend = gainindend-1;
           break;
       else
           gainindend = gainindend+1;
       end
   end
   gain(i) = str2double(tekst(gainindst:gainindend));

%    %I'm commenting this out because it creates problems for some reason
%    in some sessions. Plus it doesn't do something important apart from
%    making the script chatty
%    bkanal=tekst(modeind(1)-4:modeind(1)-1);
%    refnr=str2double(bkanal);
%    
%    if daqversion==1
%        chb(i) = 1+str2double(tekst(refindeks(refnr+1)+6:refindeks(refnr+1)+8));
%    else chb(i)=refnr+1;
%    end
   j=j+2;
end

for j=1:length(eegch)
    if eegmode(j) == 0; modeinfo='Mode A'; 
    elseif eegmode(j) == 1; modeinfo='Mode B';
    elseif eegmode(j) == 2; modeinfo='Mode -A'; 
    elseif eegmode(j) == 3; modeinfo='Mode -B';
    elseif eegmode(j) == 4; modeinfo='Mode A-B';
    elseif eegmode(j) == 5; modeinfo='Mode B-A';
    else disp('Recorded in unknown EEG mode.');
    end

    disp(sprintf('%s%i%s%i%s%s%s%i','EEG-',j,' channel ',eegch(j), '. ', modeinfo,...
        '(recorded by channel ',eegch(j),', gain ', gain(i), ').'));
end
fclose(fid2);

end


%%  SessionDuration
%
% This code is modified from EegMode.m and read session length in [sec]
% from '.set' file.
% Input (datafile) is the experiment name with the session and the extention.

%
% Takuma Kitanishi 110112

function duration = SessionDuration(datafile)

% get datafile info
daqversion=0; %old system, 1= new system
fid2=fopen(datafile,'r');
tekst=fscanf(fid2, '%c');

refindeks=strfind(tekst,'ref_');
if length(refindeks)>1
    disp('Recording was made with new system');
    daqversion=1; % new system
end

% Search 'duration' in the .set file
DurationInd = strfind(tekst,'duration ');
DurationIndSt = DurationInd + 9;
DurationIndEnd = DurationIndSt + 1;
while 1 % search end ind
    if tekst(DurationIndEnd) == sprintf('\r')
        DurationIndEnd = DurationIndEnd - 1;
        break;
    else
        DurationIndEnd = DurationIndEnd + 1;
    end
end

% Read 'duration' value
duration = str2num(tekst(DurationIndSt:DurationIndEnd));
end


%% EFG data processing functions for Fliss' part of the script

function baseline_normalisation=baseline_norm(Velocity_cmpersec_all, Fs, Fs_all, sessions, EGF_all_uV)
% DEFINING THE DATA THAT WILL BE USED FOR BASELINE NORMALISATION

%Define your baseline as the periods of immobility across the session

%index the points where the velocity is below 3cm/s

% NOTE: this removes the last 12 data points from each session
% This correponds to the last 0.48s of data.
%since this part of the code is used for calculating a baseline that will
%be used in normalisation later, the loss of half a second of data over 10
%minutes of recording is not a major concern.  This adjustment is necessary
%because of the conversions of the indexing between velocity measurements and
%frequency-power spectra. The adjusted vector is 15000 points in length.
baseline_all={};
figure
for ses=1:sessions
index_baseline=find(Velocity_cmpersec_all(:,ses)<3);    
%create a baseline index matrix, with a conversion factor of 192
baseline_index=NaN(192, length(index_baseline));
%scale the indexes from the velocity matrix so that they can be applied
%to the frequency-power matrices
for i=1:length(index_baseline)
    baseline_index(:,i)=(((index_baseline(i)*192)-191):1:(index_baseline(i)*192));
end
%reshape the baseline index into a single vector, removing duplicate
%indices
baseline_index=unique(reshape(baseline_index,[size(baseline_index,1)*size(baseline_index,2),1]));

%Create a vector the length of the EGF data
baseline_calc=zeros(length(EGF_all_uV{ses}),1);
%fill in the vector with the appropriate EGF values, as indexed by the
%baseline index
for i=1:length(baseline_index)
    baseline_calc(baseline_index(i))=EGF_all_uV{ses}(baseline_index(i));
end

%find the places where the baseline index changes by more than 1
%i.e. the prolonged periods of immobility
streak_index=find(diff(baseline_index)~=1);  
%measure the lengths of these streaks of immobility
length_streak=diff(streak_index);
length_streak=[streak_index(1); length_streak];
%parse out the periods when the animal has been still for the longest times
%the metric for this is defined as everything more than 3 standard deviations from the
%mean streak length
[immobility,immobility_index]=findpeaks(length_streak,'MinPeakHeight', (mean(length_streak)+3*std(length_streak)));

%include a quality control if statement:
% if there is only one streak of immobility, then lower the threshold for
% the length of a long immobility streak

% if you dont manage to find at least two peaks in the length of
% immobility, i.e two long streaks of the animal being still,
% then run the immobility index code again, with the threshold reduced to
% two standard deviations. 
if length(immobility_index)<2
    [immobility,immobility_index]=findpeaks(length_streak,'MinPeakHeight', (mean(length_streak)+2*std(length_streak)));
end

% find the position of the peaks by indexing the streaks by the periods of prolonged immobility 
immobility_peak=streak_index(immobility_index);
%the start of these will be the point in the streak that happens just
%before
immobility_start=immobility_index-1;
immobility_end=immobility_index;
%the data is contained within these immobility indicies
still_streak=streak_index([immobility_start,immobility_end]);

% Plot a figure to illustrate this process
subplot(2,3,ses);
%the x axis is defined as the time over all the experiments, which can be found from the sampling frequency
baseline_time=(linspace(1,length(baseline_calc)/Fs,length(baseline_calc)))';
%plot all the periods of immobiliity as our foundations of the baseline
plot(baseline_time,baseline_calc)
hold on

%plot only the prolonged periods of immobility, the top 2% longest periods,
%defined from the normal distribution
baseline_calc(1:baseline_index(still_streak(1,1)))=0;
for i=1:length(still_streak)-1
baseline_calc(baseline_index(still_streak(i,2)):baseline_index(still_streak(i+1,1)))=0;
end
baseline_calc(baseline_index(still_streak(end,2):end))=0;

%collect the prolonged immobility baseline data into a cell structure
baseline_all{ses}=baseline_calc;

plot(baseline_time, baseline_calc,'r--')
%label the axis and give the figure legends. 
ylabel('Amplitude of EEG (microvolts)'); xlabel('Total experiment time (s)');
xlim([0 max(baseline_time)]); ylim([-500 500]);
gca; set(gca,'box','off')
title(['Session ' mat2str(ses)])
end
suptitle('Indexing the baseline data for normalisation by periods of prolonged immobility');
legend('All periods of immobility', 'Prolonged periods of immobility')
legend boxoff



% we now have the relevant EFG indexed out by velocity, to include only the
% periods when the animal has a velocity of less than 3cm/s for prolonged
% periods of time
% note - 'prolonged' in this case is a relative term, meaning more than 3
% standard deviations longer than the mean period of immobility. 

% With the baseline periods defined, we can calcuate the baseline EFG.

% CALCULATING THE BASELINE 
% Here we process the baseline through the same function as is used to create the spectrogram.
% This is necessary to create comparable matrices between data and baseline
% in order for the data to be normalised against the baseline. 
BASELINE_ALL={};
for ses=1:sessions
t = (0:length(EGF_all_uV{ses})-1)'/Fs_all(ses); %I need to create time-voltage pairs for this function.
[BASELINE, logTransformed, dbconverted, t_spect, fspect] = MTSpectrogram([t baseline_all{ses}],baseline_all{ses},'show','off','frequency',Fs_all(ses), 'range',[0 100], 'tapers', [1 1]); %I used the lowest and the highest frequency I'm interested in this.
%take the mean of the baseline, across all timepoints
for j=1:size(BASELINE,2)
    if BASELINE(1,j)==0
        BASELINE(:,j)=NaN;
    end
end
% average over the processed baselines from that session to give a single
% vector for that session
BASELINE_ALL{ses}=nanmean(BASELINE,2);
end

%You now have a cell array of baseline vectors.  These can be used to normalise the data
%from each individual session, or averaged before being used to normalise
%data from all sessions. 

%The advice is to compute an average baseline, and normalise data from all
%sessions to this average. This both improves the signal to noise ratio,
%and retains he session-to-session differences. 

baseline_normalisation=mean(cell2mat(BASELINE_ALL),2);
end


function minute_by_minute_spectra(sessions, f_spect, Spectro_all)
for ses=1:sessions
    spectrogram=Spectro_all{1,ses};
    % Plot a figure of the mean power-frequency spectrum 
    % over the ten minutes of the session
    % please note: the data is smoothed via the tapering done previously
    % so the spectrograms produced are the approximate, smoothed values.
    
    %also note: the frequency goes from 4-100Hz
    % and the last 'minute' is actually 57.5 seconds 
    figure ('Name', 'Power spectra minute by minute', 'Position', [200 200 1400 700])
    hold on
    for j=1:9
        subplot(2,5,j);
        plot(f_spect(28:end),mean(spectrogram(28:end,(j*24)-23:(j*24)),2));
        ylim([0 2000])
        xlim([0 100]);
        
        set(gca, 'visible', 'off','box', 'off','xtick', []);
        legend(['Minute ' mat2str(j)]); legend boxoff
    end
    subplot(2,5,10)
    plot(f_spect(28:end),mean(spectrogram(28:end,(10*24)-23:end),2));
    ylim([0 2000])
    xlim([0 100]);
    set(gca,'visible', 'off', 'box', 'off','xtick', []);
    legend(['Minute ' mat2str(10)]); legend boxoff
    suptitle(['Session ' mat2str(ses)])
end
end


function boundedline_powerspectra(sessions, Spectro_all,f_spect)
figure ('Name', 'Power spectra', 'Position', [200 200 1400 700])
ax=zeros(sessions,1);
for ses=1:sessions
    spectrogram=Spectro_all{1,ses};
    x=(f_spect(1,28:end))';
    y=mean(spectrogram,2);
    y=y(28:end);
    e=std(spectrogram,0,2);
    e=e(28:end);
    
    % Plot a figure of the mean power-frequency spectrum
    % over the entirity of the session
    % with the standard deviation using bounded line
    % note: the frequency goes from 4-100Hz

    ax(ses)=subplot(2,3,ses);
    boundedline(x, y, e, 'b-', 'alpha');
    %ylim([0 2000])
    xlim([0 100]);
    ylabel('Power (A.U.)')
    xlabel('Frequency (Hz)');
    title(['Session ' mat2str(ses)])
    
    
end
linkaxes(ax,'xy')

end


function boundedline_normalised(sessions, Spectro_all_norm,Spectro_all,f_spect, normalisation_factor,experiment)

figure ('Name', 'Power spectra', 'Position', [200 200 1400 700])
ax=zeros(sessions,1);
for ses=1:sessions
    spectrogram=Spectro_all_norm{1,ses};
    x=(f_spect(1,28:end))';
    y=spectrogram(28:end);
    e=std(Spectro_all{1,ses},0,2);
    e=e(28:end)./normalisation_factor;
    
    % Plot a figure of the mean power-frequency spectrum
    % over the entirity of the session
    % with the standard deviation using bounded line
    % note: the frequency goes from 4-100Hz

    ax(ses)=subplot(2,3,ses);
    boundedline(x, y, e, 'b-', 'alpha');
    %ylim([0 2000])
    xlim([0 100]);
    ylabel('Power (Normalised)')
    xlabel('Frequency (Hz)');
    title(['Session ' mat2str(ses)])
    
    
end


linkaxes(ax,'xy')
ylim([0 1.5])
suptitle (['Experiment ' experiment ', normalised power spectra'])
end

%% MT spectrogram fliss Edit

%MTSpectrogram - Compute LFP spectrogram by multi-taper estimation.
%
%  The spectrogram is computed using the <a href="http://www.chronux.org">chronux</a> toolbox.
%
%  USAGE
%
%    [spectrogram,t,f] = MTSpectrogram(lfp,<options>)
%
%    lfp            wide-band LFP (one channel).
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'frequency'   sampling rate (in Hz) (default = 1250Hz)
%     'range'       frequency range (in Hz) (default = all)
%     'window'      duration (in s) of the time window (default = 5)
%     'overlap'     overlap between successive windows (default = window/2)
%     'step'        step between successive windows (default = window/2)
%     'tapers'      relative resolution and order of the tapers [NW K]
%                   (default = [3 5])
%     'pad'         FFT padding (see help for <a href="matlab:help mtspecgramc">mtspecgramc</a>) (default = 0)
%     'show'        plot spectrogram (default = 'off')
%     'parent'      parent figure or uipanel handle (default = gcf)
%     'cutoffs'     cutoff values for color plot (default = [0 13])
%    =========================================================================
%
%  NOTES
%
%    The LFP can be provided either as a time stamped matrix (list of time-voltage
%    pairs), or as a voltage vector - in which case the frequency must be specified.
%
%    The time displacement between successive short time spectra can be supplied
%    either as a 'step' (explicit time difference) or as an 'overlap' (between
%    successive time windows).
%
%  OUTPUT
%
%    spectrogram    time-frequency matrix
%    t              time bins
%    f              frequency bins
%
%  DEPENDENCIES
%
%    This function requires the <a href="http://www.chronux.org">chronux</a> toolbox.
%
%  SEE
%
%    See also CMBNSpectrogram, SpectrogramBands, PlotColorMap.

% Copyright (C) 2004-2010 by MichaÃ«l Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [spectrogram, logTransformed, dbconverted,t_spect,f_spect] = MTSpectrogram(lfp,baseline,varargin)

% Make sure chronux is installed and functional
%CheckChronux;

% Defaults
f = 1250;
frequency = [];
window = 5;
range = [];
overlap = [];
step = [];
show = 'off';
cutoffs = [0 12];
tapers = [3 5];
pad = 0;
parent = [];


% Check dependencies
if isempty(which('mtspecgramc')),
	error('This function requires the <a href="http://www.chronux.org">chronux</a> toolbox by P. Mitra, which does not appear to be installed on this system.');
end

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
end

% Check parameter sizes
if size(lfp,2) ~= 1 && size(lfp,2) ~= 2,
	error('Parameter ''lfp'' is not a vector or a Nx2 matrix (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'frequency',
			frequency = varargin{i+1};
			if ~isdscalar(frequency,'>0'),
				error('Incorrect value for property ''frequency'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'range',
			range = varargin{i+1};
			if ~isdvector(range,'#2'),
				error('Incorrect value for property ''range'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'window',
			window = varargin{i+1};
			if ~isdscalar(window,'>0'),
				error('Incorrect value for property ''window'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'overlap',
			overlap = varargin{i+1};
			if ~isdscalar(overlap,'>0'),
				error('Incorrect value for property ''overlap'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'step',
			step = varargin{i+1};
			if ~isdscalar(step,'>0'),
				error('Incorrect value for property ''step'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'tapers',
			tapers = varargin{i+1};
% 			if ~isivector(tapers,'#2','>0'),
% 				error('Incorrect value for property ''tapers'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
% 			end
		case 'pad',
			pad = varargin{i+1};
			if ~isdscalar(pad,'>-1'),
				error('Incorrect value for property ''pad'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isstring2(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
			end
		case 'cutoffs',
			cutoffs = varargin{i+1};
			if ~isdvector(cutoffs,'#2','>=0'),
				error('Incorrect value for property ''cutoffs'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'parent',
			parent = varargin{i+1};
			if ~ishandle(parent),
				error('Incorrect value for property ''parent'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).']);
	end
end

% Determine LFP frequency if it is not already provided
if isempty(frequency),
	if size(lfp,2) == 2,
		frequency = round(1/median(diff(lfp(:,1))));
	else
		frequency = f;
	end
end

% Determine step/overlap
if isempty(step),
	if isempty(overlap),
		overlap = window/2;
	end
else
	if isempty(overlap),
		overlap = window-step;
	else
		if overlap ~= window-step,
			error('Incompatible ''step'' and ''overlap'' parameters (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
		end
	end
end

% Compute and plot spectrogram
parameters.Fs = frequency;
if ~isempty(range), parameters.fpass = range; end
parameters.tapers = tapers;
parameters.pad = pad;
[spectrogram,t,f] = mtspecgramc(lfp(:,2),[window window-overlap],parameters);
t = t+lfp(1,1);
spectrogram = spectrogram';

%perform log transformation on the data
if strcmp(lower(show),'on'),
    if isempty(parent), parent = figure; end
    logTransformed = log(spectrogram);
    PlotColorMap(logTransformed,1,'x',t,'y',f,'cutoffs',cutoffs);
    colormap jet
    title ('log Transformed data')
    xlabel('Time (s)');ylabel('Frequency (Hz)');
    
    
    % dB-correction
    dbconverted = 10*log10( bsxfun(@rdivide,spectrogram,baseline));
    dbconverted(330:355,:)=repmat(min(dbconverted),26,1);                       % take out the electrical noise band
    figure
    PlotColorMap(dbconverted,1,'x',t,'y',f,'cutoffs',cutoffs);
    colormap jet
    xlabel('Time (s)');ylabel('Frequency (Hz)');
    title('dB-corrected data')
else
    logTransformed=[];
    dbconverted=[];
end
t_spect=t;
f_spect=f;
end

%% Average spectrogram function called at line 317 of the script
function[spectrogram_output]=average_spectrogram(Spectro_all,sessions, lgtrfm_all,db_all,Spectro_all_norm,normalisation_factor)
%preallocate matricies for the mean and standard deviation of the outputs
%these can be used in the calculate and plotting of average spectra with
%error bars or bounded lines. You will need the boundedline funciton to
%give bounded lines
mean_spectro=nan(length(Spectro_all{1,1}),sessions);
std_spectro=nan(length(Spectro_all{1,1}),sessions);
mean_log_transformed=nan(length(Spectro_all{1,1}),sessions);
std_log_transformed=nan(length(Spectro_all{1,1}),sessions);
mean_db_converted=nan(length(Spectro_all{1,1}),sessions);
std_db_converted=nan(length(Spectro_all{1,1}),sessions);
mean_spectro_norm=nan(length(Spectro_all{1,1}),sessions);
std_spectro_norm=nan(length(Spectro_all{1,1}),sessions);

%calculate the mean and standard deviation, session by session. 
%use the value to populate the preallocated matricies.
%perform the averaging over the columns (second dimension)
for ses=1:sessions
mean_spectro(:,ses)=        mean(Spectro_all{ses},2);
std_spectro(:,ses)=         std(Spectro_all{ses},0,2);
mean_log_transformed(:,ses)=mean(lgtrfm_all{ses},2);
std_log_transformed(:,ses)= std(lgtrfm_all{ses},0,2);
mean_db_converted(:,ses)=   mean(db_all{ses},2);
std_db_converted(:,ses)=    std(db_all{ses},0, 2);
mean_spectro_norm(:,ses)=   Spectro_all_norm{ses};                        %this value has already been averaged
std_spectro_norm(:,ses)=    (std(Spectro_all{ses},0,2))./normalisation_factor;%the standard deviation needs to be calculated from the raw data, then normalised
end

%Create an output structure
%the first branch of the structure is the type of data processing
%the second branch is the mean and standard deviation,
%organised into the 6 sessions
spectrogram_output=[];

spectrogram_output.spectro.mean=mean_spectro;
spectrogram_output.spectro.st_dev=std_spectro;

spectrogram_output.log_trfm.mean=mean_log_transformed;
spectrogram_output.log_trfm.st_dev=std_log_transformed;

spectrogram_output.db_conv.mean=mean_db_converted;
spectrogram_output.db_conv.st_dev=std_db_converted;

spectrogram_output.spectro_norm.mean=mean_spectro_norm;
spectrogram_output.spectro_norm.st_dev=std_spectro_norm;
end

%% This function is called at line 322
function [frequency_bands, theta_spect, slow_gamma_spect, fast_gamma_spect]=frequency_bands_spectra(f_spect, Spectro_all, sessions, SesDurMin)

%function to find the theta (6-11Hz), slow gamma (20-48 Hz) and fast gamma
%(70-98 Hz) data from each session, divided into minute-by minute chunks
%(10 minutes total).

%returns an output structure of frequency bands, divided by session to give
%a structure of frequency bands , and mean values for each.
% in each sub-structure, the ten columns show the 10 minutes of the
% experiment.

%PLEASE NOTE This function is originally made for the raw spectrogram data
%from 'Spectro_all', but can be used for other versions of the data, i.e.
%the decibel converted or log transformed data.  You just need to change
%the input and give the output a unique name each time the function is
%called

%make a theta index and define a theta spectrum
theta_index=(f_spect>5.9 & f_spect<11.1);
theta_spect=f_spect(theta_index);

% make a slow gamma index and define a slow gamma spectrum
slow_gamma_index=(f_spect>19.9 & f_spect<48.1);
slow_gamma_spect=f_spect(slow_gamma_index);

% make a fast gamma index and define a fast gamma spectrum
fast_gamma_index=(f_spect>69.9 & f_spect<98.1);
fast_gamma_spect=f_spect(fast_gamma_index);

%define an output cell matrix.  
% This will contain output strutures for each of the 6 sessions
frequency_bands={};

%go through the data, session by session
for ses=1:sessions
    spectro_temp=Spectro_all{ses};
    spectro_theta=spectro_temp(theta_index,:);
    spectro_slow_gamma=spectro_temp(slow_gamma_index,:);
    spectro_fast_gamma=spectro_temp(fast_gamma_index,:);
    
    
    %divide the data up minute by minute, 
    %find the average over each minute
    %do this for each of the three frequency bands
    theta_spect_minute=nan(length(theta_spect), SesDurMin);
    slow_gamma_spect_minute=nan(length(slow_gamma_spect), SesDurMin);
   fast_gamma_spect_minute=nan(length(fast_gamma_spect), SesDurMin);
     for j=1:9
        theta_spect_minute(:,j)=mean(spectro_theta(:,(j*24)-23:(j*24)),2);
         slow_gamma_spect_minute(:,j)=mean(spectro_slow_gamma(:,(j*24)-23:(j*24)),2);
         fast_gamma_spect_minute(:,j)=mean(spectro_fast_gamma(:,(j*24)-23:(j*24)),2);
    end
        theta_spect_minute(:,10)=mean(spectro_theta(:,(10*24)-23:end),2);
        slow_gamma_spect_minute(:,10)=mean(spectro_slow_gamma(:,(10*24)-23:end),2);
        fast_gamma_spect_minute(:,10)=mean(spectro_fast_gamma(:,(10*24)-23:end),2);
    
     % the data now needs to be collected into an output structure
     bands=[];
     bands.theta=theta_spect_minute;
     bands.slow_gamma=slow_gamma_spect_minute;
     bands.fast_gamma=fast_gamma_spect_minute;
     %calculate the mean value of the band for each minute of the
     %experiment within the session.
     %This is averaged minute by minute, i.e. over all rows
     bands.mean_theta=mean(theta_spect_minute);
     bands.mean_slow_gamma=mean(slow_gamma_spect_minute);
     bands.mean_fast_gamma=mean(fast_gamma_spect_minute);
     
     %the structure is saved to the output cell array
     frequency_bands{ses}=bands;
end
end


%% fft function that is used in line 132


function xf = fftbandpass(x,Fs,Fs1,Fp1,Fp2,Fs2)
% function XF = fftbandpass(X,FS,FS1,FP1,FP2,FS2)
%
% Bandpass filter for the signal X (time x trials). An acuasal fft
% algorithm is applied (i.e. no phase shift). The filter functions is
% constructed from a Hamming window.
%
% Fs : sampling frequency
%
% The passbands (Fp1 Fp2) and stop bands (Fs1 Fs2) are defined as
%                 -----------
%                /           \
%               /             \
%              /               \
%             /                 \
%   ----------                   -----------------
%           Fs1  Fp1       Fp2  Fs2
%
% If no output arguments are assigned the filter function H(f) and
% impulse response are plotted.
%
% NOTE: for long data traces the filter is very slow.
%
%------------------------------------------------------------------------
% Ole Jensen, Brain Resarch Unit, Low Temperature Laboratory,
% Helsinki University of Technology, 02015 HUT, Finland,
% Report bugs to ojensen@neuro.hut.fi
%------------------------------------------------------------------------

%    Copyright (C) 2000 by Ole Jensen
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You can find a copy of the GNU General Public License
%    along with this package (4DToolbox); if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

if size(x,1) == 1
    x = x';
end
% Make x even
Norig = size(x,1);
if rem(Norig,2)
    x = [x' zeros(size(x,2),1)]';
end

% Normalize frequencies
Ns1 = Fs1/(Fs/2);
Ns2 = Fs2/(Fs/2);
Np1 = Fp1/(Fs/2);
Np2 = Fp2/(Fs/2);

% Construct the filter function H(f)
N = size(x,1);
Nh = N/2;

B = fir2(N-1,[0 Ns1 Np1 Np2 Ns2 1],[0 0 1 1 0 0]);
H = abs(fft(B));  % Make zero-phase filter function
IPR = real(ifft(H));
if nargout == 0
    subplot(2,1,1)
    f = Fs*(0:Nh-1)/(N);
    plot(f,H(1:Nh));
    xlim([0 2*Fs2])
    ylim([0 1]);
    title('Filter function H(f)')
    xlabel('Frequency (Hz)')
    subplot(2,1,2)
    plot((1:Nh)/Fs,IPR(1:Nh))
    xlim([0 2/Fp1])
    xlabel('Time (sec)')
    ylim([min(IPR) max(IPR)])
    title('Impulse response')
end


if size(x,2) > 1
    for k=1:size(x,2)
        xf(:,k) = real(ifft(fft(x(:,k)) .* H'));
    end
    xf = xf(1:Norig,:);
else
    xf = real(ifft(fft(x') .* H));
    xf = xf(1:Norig);
end

end


%% Extract oscilation phase from the band pass filtered EEG signal by Hilbert transform. Troughs are set to 0/360 deg. Output is in deg
function [Phase, Amp] = Osciphase(Eeg)

% Hilbert transform
Z = hilbert(Eeg);

% Wave amplitude
Amp = abs(Z);

% EEG phase in rad
Phase = angle(Z);

% Rad to Degree (-180 to +180)
Phase = Phase / pi *180;

% Degree (0 to +360)
Phase = Phase + 180;

end

%% Detect gamma events and returns peak index
% Power, EEG power vector
% Fs, EEG sampling frequency [Hz]
% wLength, gamma temporal window around peaks [sec]
% peakTiming, Peak positions thresholded with mean+2SD [sec]
% gammaLength, Temporal length within gamma events [sec]
function [peakTiming, gammaLength] = gammaEvent(Power,Fs,wLength)

% Threshold (mean + 2 SD)
Th = mean(Power) + 2 .* std(Power);

% Local maximum
[pks,locs] = findpeaks(Power);

% Select peak index
locs = locs(pks>Th);

% Remove initial and last index
if locs(end) == length(Power)
    locs = locs(1:end-1);
end
if locs(1) == 1
    locs = locs(2:end);
end

% Conversion from eeg index to timing [sec]
peakTiming = locs./Fs;


% Half-window length [in idx unit]
wLength = wLength./2;
wLength = round(wLength.*Fs);

% Initialitze tIdx (=index in gamma events)
tIdx = zeros(length(Power)+2*wLength,1);

% Search index in gamma events
for i=1:length(locs)
    idx = (locs(i)-wLength:locs(i)+wLength)+wLength;
    tIdx(idx)=1;
end
tIdx = tIdx(wLength+1:length(Power)+wLength);

% Temporal length in gamma events [sec]
gammaLength = sum(tIdx)./Fs;
end


%% Finds the thesholds for 1st, 2nd, 3rd and 4th interquartile (IQ) interval

%First quartile: the lowest 25% of numbers <q1
%Second quartile: between 25.1% and 50% (up to the median) q1<
%>q2(essentially median)
%Third quartile: 51% to 75% (above the median) q2< >q3
%Fourth quartile: the highest 25% of numbers >q3

function [q1, q2, q3]=findiquartile(V)


% Check data is columnwise
if size(V,1)<size(V,2)
    V = V';
end
    
V_sort=sort(V,'descend'); %sort the values in the input vector to get the top range

q1 = max(V_sort(ceil(length(V)*0.75):end));

q2 = median(V_sort);

q3 = min(V_sort(1:ceil(length(V)*0.25)));
end

%% Select spikes that happend in Timing[sec]+-wLength/2[sec]
% ts, spike timing[sec]
% Timing, event timing such as gamma event peaks [sec]
% wLength, full window length [sec]
% idx, selected spike index
function idx = selectSpikes(ts, Timing, wLength)

% Spike number
N = length(ts);

% Initialization
idx = zeros(N,1);

% Half-window length
wLength = wLength./2;

% For each spike
for i=1:N
    if any(ts(i)-wLength<Timing & Timing<ts(i)+wLength)
        idx(i)=1;
    end
end
end