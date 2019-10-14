%% Written by Antonis Asiminas, SIDB July 2019.
% Updated September 2019
% This script makes use of functions from FMAT, chronux, PAC (Onslow et al 2010) and the circular statistics toolbox as well as some functions used in/written for Kitanishi et al 2015 Neuron and Tort et al 2008 PNAS 

% The inputs are: 
% mypos position data
% .set file from Axona system
% .egf files recorded in Axona system which contain 4800Hz continuous LFP data, 
% an xl file containing cluster identities from previously analysed data,
% .clu and .fet files from Klusters were I get the cluster spike times.
% These are from previsouly combined sessions in order to spike sort the same clusters accross the experiemntal paradigm.  

% The clusters used in this analysis have been chosen because they are
% putative pyramidal cells nicely isolated from noise. They are the same
% clusters that have been use in SpikeTrainAnt

% The way I run this is pretty silly. I run it separatelly for every
%rat/session and then get the outputs, average across genotypes and run
%stats in SPSS or R and create graphs in Graphpad/Inkscape. In the future I
%will need to write a master script that is erading a .txt file with the
%absolute directories of all the sessions I'm interested in. 

%I have tried to comment as much as possible throughout the script. I have no doubt there will be bugs, redundant variables, optimization issues etc. 
%The way this code is written leads to a very messy workspace in the end
%which makes debugging challenging.
%If you have any questions please get in touch in a.asiminas@ed.ac.uk and
%I'll try to help ask much as I can. Remember this is the second analysis
%script I've ever written, be understanding.

% This script is separated into five sections:

% In the fist section the position data are imported, separated in sessions and put into
% usable matrices

% In the second section the LFP data are imported and analysed. Basic
% spectrograms and the power of oscillatory bands of interest are
% calculated different subsets of the data (timewise and velocity-wise).
% Moreover certain metrics necesary for subsequent sections are being
% computed

% In section three I focus on Phase-amplitude coherence. This is computed for entire sessions and on a minute by minute and its recurrence
% over the experiment is explored by correlating the comodulogram matrices.

% In the fourth section the spike timestamps are imported and are corrected
% if they're coming from previously combined sessions (session number >1)

% In the fith section the LFP and spike timestamps are coming together in
% order to calculate phase locking of spikes to different oscillatory
% bands (theta, slow gamma and medium/fast gamma). 

% The last section is about taking everything together and creating plots and output tables for further analysis.


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
%Dvorak et al 2018 PLoS Biol and Colgin et al 2009 Nature. 
% I could have done that with a struct but it's difficult to index if I need to
% Theta bandpass
Lims{1,1} = 5; %low stop for theta
Lims{1,2} = 6; %low pass for theta
Lims{2,1} = 11; %high pass for theta
Lims{2,2} = 12; %high stop for theta
% Slow gamma bandpass
Lims{1,3} = 26; %low stop for slow gamma
Lims{1,4} = 27; %low pass for slow gamma
Lims{2,3} = 48; %high pass for slow gamma
Lims{2,4} = 50; %high stop for slow gamma
% Fast gamma bandpass
Lims{1,5} = 58; %low stop for fast gamma
Lims{1,6} = 60; %low pass for fast gamma
Lims{2,5} = 98; %high pass for fast gamma
Lims{2,6} = 100; %high stop for fast gamma

% Gamma event window length [sec]
% Spikes which happend in (gamma power peak)+-(wLength/2) will be used.
% Default, 0.400 sec according to Colgin et al., 2009 Nature.
wLength = 0.400;
pixelratio=330; %that's pixels per meter. This has to be determined manually by the experimenter.

startup; %this calls a function at the bottom of the script for figure uniformity

%% Section One %%%% Loading position data that will be used for excluding spikes and LFP during imobility.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Reading positions ...');

pos_all=load('merge.mypos'); %I'm loading the merged position file. This is better that the individual session mypos files bcz it has gone through Kalman filtering.

%Extract series of coordinates
pos_x_all = pos_all(:,1);
pos_y_all = pos_all(:, 2);
pos_time_all = pos_all(:,4)*100000; %I multiply in order to be able to use the Tseg file to split
disp('...done');

T_segs =load('merge.tseg'); %use the tseg file to see the size of the segments and the number of them
T_segs=T_segs(:,2);
T_segs = [T_segs(1:end);T_segs(end)+(T_segs(end)-T_segs(end-1))]*100000; %that is to correct with the spikes later

%this separates the x-y-t data into sessions
disp('Analysing movement...');
Pos_X_all_ses = zeros(SesDurMin*60*50, sessions); %that's based on the assumption that all session should be the same length. Probably not the most elegant way to write that. If they were not the same I should use a cell array anyway
Pos_Y_all_ses = zeros(SesDurMin*60*50, sessions);
Pos_Time_all_ses = zeros(SesDurMin*60*50, sessions); %I need 600sec*50Hz = 30000 samples for each of the sessions. Probably not the most elegant way to write that either.

for tp=1:length(T_segs)-1 
    vals=[]; %clean vals after every iteration
    ind_start=find(pos_time_all==T_segs(tp));
    ind_dif=find(pos_time_all==T_segs(tp+1))-find(pos_time_all==T_segs(tp));
    if tp==6
        ind_dif=length(pos_time_all)-find(pos_time_all==T_segs(tp));
    end
    
    if ind_dif>SesDurMin*60*50  %usually a 10min session in axona is about 30049 datapoints
        n=1;
        while n<SesDurMin*60*50+1    %here I take only the datapointsthat correspond to 50Hz*Session in seconds
        vals(n,1) = pos_time_all(ind_start);
        ind_start=ind_start+1;
        n=n+1;
        end 
        
    Pos_Time_all_ses(:,tp)=vals-T_segs(tp); %If everything is right then all columns of Pos_Time_all_ses should be the same. I also realised that the system has a lag of about 1sec (0.93sec). That means that for that period the EEG is 0 and nothing is happening. I want correct that by taking only the session time. So I need 600sec*50Hz = 30000 samples for each of these   
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

%% Section Two %%%%%%%%%%%%%%%%% Reading EGF file and doing some basic analyses of the LFP data. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This section will provide a number of outputs suitable for plotting but
% also some variables that will be used together with the spiketiming data
% (Section Four)for the phase locking analysis.



%Reading EGF files and their sampling frequencies from the folder and
%adding them in cell arrays
EGF_all = cell(1,sessions); %a cell array with the EEG for all sessions
Fs_all = zeros(1,sessions); %a vector with the sampling rate of all sessions' EEG
for ses=1:sessions
    egffile = [experiment, num2str(ses),'.egf'];
    disp(['Reading eeg data from session ', num2str(ses)])
    [EGF,Fs] = readEGF(egffile); %function at the bottom of the script
    if length(EGF)< Fs*60*SesDurMin+ 4800   %The eeg files are supposed to have 1 more second worth of samples that are 0. I'just try to have uniform data and then I'll remove it 
        
        while length(EGF)-(Fs*60*SesDurMin + 4800)<0 %for my analysis later I only take the samples for Fs*60*SesDurMin (for a 10min session in 4.8kHz sampling rate it's 2880000)
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
                disp(['Recorded in -Ref mode! EEG sign was corrected in session ', num2str(ses)])
    end
    eegmode_all(ses)=eegmode;
    Gain_all(ses)=Gain;
end


% Convert EEG values to micro Volt (uV) see the email from Jim below
EGF_all_uV=cell(1,sessions); %a cell array with the EEG for all sessions in micro Volt values
for ses=1:sessions
EGF_all_uV{ses} = EGF_all{ses}(1:Fs*60*SesDurMin) ./ 32768 .* 1500000 ./ Gain_all(ses); % these are stored at 16-bit resolution, so the values range from -32768 to 32767
end

EGF_all_uV_Ztransformed=cell(1,sessions); %with this I z-tranform my EEG rwa signal. I'll calculate the power spectrum densities based on that as well. 
for ses=1:sessions
    EGF_all_uV_Ztransformed{ses}=(EGF_all_uV{ses}-mean(EGF_all_uV{ses}))./std(EGF_all_uV{ses});
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


%% This is to reshape the long EEG signal of every session to 500ms non-overlapping epochs. I also take only the first 600sec for each session and not the lag 0.93s at the end
% From these segments I can export some examples. I just need to run the fft bandpass filters on a good example section(s)

EGF_500ms_uV_all=cell(1, sessions);
EGF_500ms_uV_Ztransformed_all=cell(1, sessions);
for ses=1:sessions
EGF_500ms=reshape(EGF_all_uV{ses}(1:Fs*60*SesDurMin), [Fs*0.5 60*SesDurMin/0.5]);
EGF_500ms_Ztransformed=reshape(EGF_all_uV_Ztransformed{ses}(1:Fs*60*SesDurMin), [Fs*0.5 60*SesDurMin/0.5]);
EGF_500ms_uV_all{ses}=EGF_500ms;
EGF_500ms_uV_Ztransformed_all{ses}=EGF_500ms_Ztransformed;
end

%This is time in seconds. Time_500ms_all has all the sampling timepoints for EEG in a given session 
Time_ses=(0:Fs*60*SesDurMin-1)'/Fs +1/Fs; %that also corrects for the time of the first sampling with should be 1/Fs
Time_500ms_all=reshape(Time_ses,[Fs*0.5 1200]); %this will be used to index time segments. =< from first row and >from last row of each segment (column). The element in between could be thought as redundant but 

% Calculate speed for all the epochs. I can export that for the movement analysis to show that genotypes did not different in speed of movement.
% Each row is a session and it's column is the velocity for a given 500ms
% epoch that corresponds to the different columns in every cell of EGF_500ms_uV_all
Velocity_cmpersec_500=NaN(sessions,60*SesDurMin/0.5);  
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


% I identify here all the periods that the animal is moving above the
% threshold for 2s
Velocity_threshold = 3; %that's cm/s and I can change that if I need different threshold
EGF_uV_2s_movement=cell(1, sessions);
EGF_uV_ztransf_2s_movement=cell(1, sessions);
for ses=1:sessions
disp(['Finding periods of immobility and excluding EEG in session ', num2str(ses)]);
    % Logical vector where velocity > threshold. I'll use that to take the
    % bits from the EEG signal
    binaryVector = Velocity_cmpersec_500(ses,:) > Velocity_threshold;
    % Label each region with a label - an "ID" number. The next few
    % functions are from the Image processing toolbox in Matlab
    [labeledVector, numRegions] = bwlabel(binaryVector);

    % Measure lengths of each region and the indexes. The output hear is a
    % structure that has all the bits from the vector that are above the
    % threshold, organised based on the concecutive elements.
    measurements = regionprops(labeledVector, Velocity_cmpersec_500(ses,:), 'Area', 'PixelValues');

    % Find regions where the area (length) are 2s or greater and
    % put the values into a cell of a cell array
    movement_epochs=cell(1,numRegions);
for k = 1 : numRegions
  if measurements(k).Area >= 4 %that is because my ephochs are 500ms therefore 2 seconds are 4xepochs
    % Store the values in movement_epochs. This eventually si have empty
    % cells and cells that have more than 4 elements in them
    movement_epochs{k} = measurements(k).PixelValues;
  end
end
% Display the regions that meet the criteria:
count=1; %start a counter to use as an index for all the 2s movement epochs I will get from each session. It'll reset to 1 when I switch to the next session
for k = 1 : numRegions % go through the movement_epochs with contains the velocities corresponding to more than 4 500ms epochs
    if ~isempty(movement_epochs{k}) %take only the cells that are not empty; in other words contain the above threshold speed for 4*500ms epochs
        %if length(movement_epochs{k})>=4 %only the cells that have more than 2s worth of movement above threshold. By definition there should not be any cells with vectors shorter than 4
        
            times=fix(length(movement_epochs{k})/4); %check if the period is multiples of 4s Of course I'm losing some data here from all the movement periods that are not multiples of 2s. But that's ok.
            
            for i=1:4:times*4 %take now all the 2s chunks for this particular movement_epoch
                %here i'm finding the index of the start of the 2s epoch. 
                index=strfind(Velocity_cmpersec_500(ses,:),movement_epochs{k}(i:i+3));
                %I now use this to index the EEG files (uV.....
                EGF2s=EGF_500ms_uV_all{ses}(1:end,index:index+3); %the structure of EGF_500ms_uV_all{ses} is that every column has the EEG data for a 500ms epoch
                %....and z-transformed)
                EGF2s_ztrans=EGF_500ms_uV_Ztransformed_all{ses}(1:end,index:index+3); %EGF_500ms_uV_Ztransformed_all has the same structure as EGF_500ms_uV_all
                
                %here I use reshape bsz I want to have each 2s chunk as a
                %column. It'll help with the spectrogram calculation
                EGF_uV_2s_movement{ses}(1:Fs*2,count)=reshape(EGF2s,[size(EGF2s,1)*size(EGF2s,2) 1]);
                EGF_uV_ztransf_2s_movement{ses}(1:Fs*2,count)=reshape(EGF2s_ztrans,[size(EGF2s_ztrans,1)*size(EGF2s_ztrans,2) 1]);
                
                count=count+1; %if that works as is should be it'll count up every time a new 2s chunk is identified 
            end
    end
end

end


disp('Calculating EEG power in entire sessions and 2min chunks without excluding periods of immobility');
%This is to caclulate the spectrogram and band powers for the entire sessions without
%excluding periods of immobility (<3cm/s). I will then chop into 2min
%epochs calculating the individual band power in these epochs
Spectrogram_sessions_all=NaN(86,12); %The first 6 columns are spectrogrmam based on uV values and the next 6 are for spectrogram based on z transformed 
Spectrogram_log_sessions_all=NaN(86,12); %Im cheating here because I know that the length of f_spect_ses is 86 for 0-100Hz. This is not elegant of course and I need to find another way for preallocation
Bands_2mins=cell(1, sessions); %every cell has the power of three bands calculated from uV and z-transformed EEG data for five 2min periods of each session. The first column is the average velocity for that 2min period
Bands_all_ses=NaN(3,12); %We look at 3 bands (Theta, slow gamma and medium gamma) for 6 sessions (from uV and z-transformed data).
Bands_500bins=cell(3,2); %Every cell has the power of a band for all 500bins for all sessions; Cell array rows are freq bands. The first array column is spectrogrm fromEEG uV and the second is from z-trnsfrm EEG
vbin_borders=0:2:40;%these are the velocity bins I use
vbin_middle=1:2:39; %that has the middle of the velocity bins so I can plot after
BandsPWR_Velbins=cell(3,2); %the rows are freq bands and the colums are one normal one z-transf EEG power. Each cell have the power of a band for all vbin velocity for all 6 sessions
for ses=1:sessions
    %I'm calculating the spectrogram for every 2s moving period of each
    %session and then average accross all of them to get the power
    %spectrogram for that session/rat. 
    %the window is 500ms wide and I have 95% overlap between windows
    [spectrogram, logTransformed,~,~] = MTSpectrogram([Time_ses EGF_all_uV{ses}],'frequency',Fs_all(ses),'window', 0.5,'overlap',0.475, 'range',[0 100]);
    
    [spectrogramZ, logTransformedZ,t_spect_ses,f_spect_ses] = MTSpectrogram([Time_ses EGF_all_uV_Ztransformed{ses}],'frequency',Fs_all(ses),'window', 0.5,'overlap',0.475, 'range',[0 100]);
    
    Spectrogram_sessions_all(:,ses)=mean(spectrogram,2);
    
    Spectrogram_sessions_all(:,ses+6)=mean(spectrogramZ,2);
    
    Spectrogram_log_sessions_all(:,ses)=mean(logTransformed,2);
    
    Spectrogram_log_sessions_all(:,ses+6)=mean(logTransformedZ,2);
    
    f_lims=[6 11;25 44;52 86]; %these are the positions in vector f_spect_ses for theta, slow gamma and fast gamma (found manually)
    
    for i=1:length(f_lims)
    
        for epoch=1:size(Velocity_cmpersec_500,2)
        
            [~,idx_start]=min(abs(t_spect_ses-Time_500ms_all(1,epoch)));
        
            [~,idx_end]=min(abs(t_spect_ses-Time_500ms_all(end,epoch)));
   
            Bands_500bins{i,1}(ses,epoch) = mean(spectrogram(f_lims(i,1):f_lims(i,2), idx_start:idx_end),'all');
            Bands_500bins{i,2}(ses,epoch) = mean(spectrogramZ(f_lims(i,1):f_lims(i,2), idx_start:idx_end),'all');
            
        end
    end
    
    
    disp(['Calculating how EEG power is modulated by speed in session ', num2str(ses)]);
    %Here I bin the velocity into 20bins and calculate the power for all
    %three freq bands of interest
    for band=1:length(f_lims)
        for v=1:length(vbin_middle)     
            %this gives the mean power for all the times the velocity was
            %between the different borders of vbin_borders. It makes
            %plotting easier
            BandsPWR_Velbins{band,1}(ses,v)=mean(Bands_500bins{band,1}(ses, Velocity_cmpersec_500(ses,:)>=vbin_borders(v) & Velocity_cmpersec_500(ses,:)<vbin_borders(v+1)));
            BandsPWR_Velbins{band,2}(ses,v)=mean(Bands_500bins{band,2}(ses, Velocity_cmpersec_500(ses,:)>=vbin_borders(v) & Velocity_cmpersec_500(ses,:)<vbin_borders(v+1)));
       
        end
    end
    

    %here I calculate the power for the bands of interest without excluding
    %any periods of low velocity (<3cm/s)
    for i=1:length(f_lims)
    
        Bands_all_ses(i,ses)= mean(spectrogram(f_lims(i,1):f_lims(i,2),:),'all');
        Bands_all_ses(i,ses+6)=mean(spectrogramZ(f_lims(i,1):f_lims(i,2),:),'all');
        
    end
    
    %I'm now splitting the spectrograms in 2min epochs to see if there are any patterns in the power of specific bands of interest (Theta, slow gamma, fast gamma).
    t_lims=[1 4791;4792 9591;9592 14391;14392 19191;19192 23981]; %these are the positions in vector t_spect-ses for 5 2min periods within each session (found manually)
    
    
    for i=1:length(f_lims) %for all the bands
        for j=1:length(t_lims) %for all the 2min periods
           
             Bands_2mins{ses}(j,i+1)= mean(spectrogram(f_lims(i,1):f_lims(i,2),t_lims(j,1):t_lims(j,2)),'all'); %that is to capture the power of all bands in every 2min period of a session
             Bands_2mins{ses}(j,i+4)= mean(spectrogramZ(f_lims(i,1):f_lims(i,2),t_lims(j,1):t_lims(j,2)),'all');
        
             [~, col_start]=find(Time_500ms_all==t_spect_ses(t_lims(j,1))); %finding with velocity 500ms bins correspond to the start and end of the 2min period
             [~, col_end]=find(Time_500ms_all==t_spect_ses(t_lims(j,2)));
             Bands_2mins{ses}(j,1)=mean(Velocity_cmpersec_500(ses,col_start:col_end),2); %tht gives me the average velocity for every 2min period. Since I have not excluded the immobility in this calculation
             
        end
    end
 
end   


disp('Calculating EEG power in entire sessions excluding periods of immobility');
% Here I calculate the spectrograms for every session only including periods
%of mobility I have Identified previously in loop in line280. The way I
%did that means that some seconds of time will be not being used but that's
%not a big problem for the current question
t = (1/Fs:1/Fs:2)'; %this is the time vector for the 2sec movement epochs. This doesn't correspond to real time since the 2s movement periods are from all over the sessions.
Spectrogram_move=NaN(length(f_spect_ses),12); %The first 6 columns are spectrogrmam based on uV values and the next 6 are for spectrogram based on z transformed 
Spectrogram_log_move=NaN(length(f_spect_ses),12);
Bands_move_ses=NaN(3,12); %We look at 3 bands (Theta, slow gamma and medium gamma) for 6 sessions (from uV and z-transformed data). 
for ses=1:sessions
    epoch_spectgram=NaN(length(f_spect_ses),size(EGF_uV_2s_movement{ses},2)); %f_spect_ses is the same as f_spect since I'm using the same frange [0 100]
    epoch_log_spectgram=NaN(length(f_spect_ses),size(EGF_uV_2s_movement{ses},2));
    epoch_Z_spectgram=NaN(length(f_spect_ses),size(EGF_uV_2s_movement{ses},2));
    epoch_Z_log_spectgram=NaN(length(f_spect_ses),size(EGF_uV_2s_movement{ses},2));
    for mov_2s_epoch=1:size(EGF_uV_2s_movement{ses},2)
        %the window here as before is 500ms wide and I have 95% overlap
        %between windows
        [spectrogram, logTransformed,~,~] = MTSpectrogram([t EGF_uV_2s_movement{ses}(:,mov_2s_epoch)],'frequency',Fs_all(ses),'window', 0.5,'overlap',0.475, 'range',[0 100]); 
        epoch_spectgram(:,mov_2s_epoch)=mean(spectrogram,2);
        epoch_log_spectgram(:,mov_2s_epoch)=mean(logTransformed,2);
        %this is calculting the exact same thing but for the
        %z-transformed data EEG data
        [spectrogramZ, logTransformedZ,t_spect,f_spect] = MTSpectrogram([t EGF_uV_ztransf_2s_movement{ses}(:,mov_2s_epoch)],'frequency',Fs_all(ses),'window', 0.5,'overlap',0.475, 'range',[0 100]); 
        epoch_Z_spectgram(:,mov_2s_epoch)=mean(spectrogramZ,2);
        epoch_Z_log_spectgram(:,mov_2s_epoch)=mean(logTransformedZ,2);  
    end
    
    Spectrogram_move(:,ses)=mean(epoch_spectgram,2);
    Spectrogram_move(:,ses+6)=mean(epoch_Z_spectgram,2);
    Spectrogram_log_move(:,ses)=mean(epoch_log_spectgram,2);
    Spectrogram_log_move(:,ses+6)=mean(epoch_Z_log_spectgram,2);
    Bands_move_ses(1,ses)=mean(Spectrogram_move(6:11,ses)); %Theta
    Bands_move_ses(1,ses+6)=mean(Spectrogram_move(6:11,ses+6));
    Bands_move_ses(2,ses)=mean(Spectrogram_move(25:44,ses)); %Slow gamma
    Bands_move_ses(2,ses+6)=mean(Spectrogram_move(25:44,ses+6));
    Bands_move_ses(3,ses)=mean(Spectrogram_move(52:end,ses)); %Medium gamma
    Bands_move_ses(3,ses+6)=mean(Spectrogram_move(52:end,ses+6));
   
end




%% 
disp('Choosing a nice EEG example for every session with strong theta, slow gamma and fast gamma');
% This is to find some nice example traces for figures. I NEED TO PLOT
% THESE EXAMPLES AS SVG AND THEN USE AN EXAMPLE FROM EACH SESSION IN FIGURE
% 3 OF THE PAPER
% I'm finding the epochs when I have nice theta and gamma power
StrongLFP_bandpass=cell(1+length(Lims)/2,sessions); %this will have 4 example traces for each session from that rat. Top is raw, second is theta, third is slow gamma and fourth is fast gamma
for ses=1:sessions
    
    [Tq1, Tq2, Tq3] = findiquartile(Bands_500bins{1,1}(ses,:)); %I only use the powers calculated from the EEG uV NOT the z-scored
    [Sq1, Sq2, Sq3] = findiquartile(Bands_500bins{2,1}(ses,:));
    [Fq1, Fq2, Fq3] = findiquartile(Bands_500bins{3,1}(ses,:));
    
    Ttop= Bands_500bins{1,1}(ses,:)>=Tq1; %depending what the findings from the band power are, I can change that a get a medium power theta or gamma to illustrate a point
    Stop= Bands_500bins{1,1}(ses,:)>=Sq1;
    Ftop= Bands_500bins{1,1}(ses,:)>=Fq1;
    
    Alltop = Ttop & Stop & Ftop; %That is a logical vector. Ones indicate the 500ms epochs where theta, slow gamma and high gamma are strong.
    
    StrongLFP=EGF_500ms_uV_all{ses}(:,Alltop);
    indx=2;
    for i=1:2:length(Lims)-1
        %I'm taking the middle segment ceil(end/2) from all the segments
        %that meet the criterion in line 424
        StrongLFP_bandpass{1,ses}=StrongLFP(:,ceil(end/2)); %getting the raw LFP signal
        StrongLFP_bandpass{indx,ses} = (fftbandpass(StrongLFP(:,ceil(end/2)),Fs, Lims{1,i}, Lims{1,i+1}, Lims{2,i}, Lims{2,i+1}))'; % Theta, slow-gamma and fast/medium-gamma bands         
        indx=indx+1;
    end
    
end



%% -Here I calculate instateneous Amp and Phase that I'll be using in the Phase Amplitude Coherence analysis (Section Three) and Spike phase locking analysis (Section Four)------------------------------------------------------------------------------%%%%%%%%%%%%%%%%%%%

% Band pass filter. I use the frequency limits set at the start of this
% script and apply them to the raw LFP for as many times as the frequency
% bands I'm interested in. In this case I'm interested in theta slow-gamma
% and medium/fast gamma. From this point on I will only use the LFP in uV
% NOT the z-transformed
%this is a cell array with 3 rows 6 columns. Each row correspond to an LFP band and each column corresponds to a session. Top row is theta, second row is slow gamma and third row is medium/fast gamma.

EGF_bandpass_all= cell(length(Lims)/2, sessions); 
for ses=1:sessions
    indx=1;
    for i=1:2:length(Lims)-1
    EGF_bandpass = fftbandpass(EGF_all_uV{ses}, Fs, Lims{1,i}, Lims{1,i+1}, Lims{2,i}, Lims{2,i+1}); % Theta, low-gamma and high-gamma bands
    EGF_bandpass_all{indx,ses}=EGF_bandpass'; %I invert just because I like column vectors
    indx=indx+1;
    end
end

% Calculate eeg instantaneous phase based on TROUGHs (that is 0 phase means though) with Hilbert transform
disp('Computing instantaneous Phase & Amplitute with Hilbert transform for all sessions and oscillation bands')
Phase_all=cell(length(Lims)/2,sessions);
Phase_rad_all=cell(length(Lims)/2,sessions); %I will need that for PAC (Section Three)
Amp_all=cell(length(Lims)/2, sessions);
for ses=1:sessions
    for i=1:length(Lims)/2
        [phase,Phase_rad,Amp] = Osciphase(EGF_bandpass_all{i,ses});
        Phase_all{i,ses}=phase;
        Phase_rad_all{i,ses}=Phase_rad;
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


% Detect peak timing & temporal length of gamma events [sec]
%the name of the function is misleading since I'm detecting events for all
%oscillatory bands of interest. Not sure where will I use the theta events though.
%The slow and fast gamma events will be used on the final stage of this script where spikes and LFP will come together.
%Kitanishi et al 2015 Neuron and Colgin et al Nature 2009 only use spike during strong gamma events  to
%calculate phase locking
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


disp('Computing Theta and Gamma Phase assymetries for every minute of the experiment')
% This is related to Trimper et 2014 Hippocampus (Fig6) and Kitanashi et al
% 2015 Neuron (Fig2H. Here I use 20deg bins (uniform distribution 5.56%). KItanishi et al used 30deg bins (uniform
% distribution 8.3%)
Osc_Ph_bins=0:20:360;
Phase_portion_all=cell(3,6); Phase_portion_all(:)={NaN(10,length(Osc_Ph_bins)-1)}; %every cell has the portions of phases for give oscillation (row) phase bin in the 10min of a given session (column)
min_borders=0:length(Phase_all{1,1})/10:length(Phase_all{1,1}); %this is to define the minute borders. It'll be usefull for the PAC as well and any other min by min measure
for ses=1:sessions
    for band=1:length(Lims)/2     
        for minut=1:length(min_borders)-1
            for ph=1:length(Osc_Ph_bins)-1
                %this calculates the porportion of a certain phase of an
                %oscillation compared to the entirety of the signal. For 90deg bin
                %that should be .25 (Trimper et al 2014) if the wave was
                %symmetrical
                Phase_portion_all{band,ses}(minut,ph)= length(Phase_all{band,ses}(Phase_all{band,ses}(min_borders(minut)+1:min_borders(minut+1))>=Osc_Ph_bins(ph) & Phase_all{band,ses}(min_borders(minut)+1:min_borders(minut+1))<Osc_Ph_bins(ph+1)))/length(Phase_all{band,ses}(min_borders(minut)+1:min_borders(minut+1)));    
            end
            
        end
    end
end


%% Section Three %%%%%%%%%%%%%%%%% Analysis of Phase-Amplitude Coupling between Theta and Gamma. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%In previous versions of this I have used functions from the PAC toolbox
%(Onslow et al 2010). I realised that this might not be the best approach
%(DON'T ask me why!) so here I'm using the approach decribed previously in 
%Tort et al 2008 PNAS.



%this is for minute by minute PAC analysis. Checking PAC in the entire session is
%meaningless because the temporal dynamics of this phenomenon/metric are
%very fast. I will average though in half sessions to have an adittional
%output

% Defining Theta phase bins

nbin = 18; % number of phase bins. Here I use 18 20deg bins like in Tort et al.
Ph_pos=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin 
    Ph_pos(j) = -pi+(j-1)*winsize; 
end
disp('Computing PAC for all Theta Power interquartile (IQ) intervals for all minutes of the experiment')


Gamma_MeanAmp=cell(2,4); %first row is for slow gamma and secod is for fast gamma. The four columns are for the 4 theta power IQs. Every cell will have all minutes
SG_ModIndex=NaN(4,60); %row are for IQs, columns are for minutes
FG_ModIndex=NaN(4,60);


countmins=1; %counter I need for organising the output
for ses=1:sessions 
    for minut=1:length(min_borders)-1 %here I'll be indexing a vector created in line 565
    
        [q1, q2, q3] = findiquartile(Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))); %finds the thresholds for 1st, 2nd, 3rd and 4th interquartile (IQ) interval of Theta power. This is based on Tort et al 2009 PNAS (Fig 4A more specifically)
        
        Phase_1IQi=Phase_rad_all{1,ses}(Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))<q1); %this gives a vector with the instantaneous phase of theta during the time theta power is at it's the weakest
        SG_Amp_1IQi=Amp_all{2,ses}(Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))<q1); %these two give vectors with the instantaneous amplitude of slow and fast gamma during the same instances as theta phases. They have the same length as Phase_1IQi.
        FG_Amp_1IQi=Amp_all{3,ses}(Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))<q1);
        
        [MI_SG_1IQi,MeanAmps_SG_1IQi]=ModIndex(Phase_1IQi, SG_Amp_1IQi, Ph_pos); %I use these to calculate the modulation index and a vector of mean gamma amplitudes for every theta phase bin (same length as nbin)
        [MI_FG_1IQi,MeanAmps_FG_1IQi]=ModIndex(Phase_1IQi, FG_Amp_1IQi, Ph_pos);
        
        Phase_2IQi=Phase_rad_all{1,ses}(Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))>=q1 & Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))<q2);
        SG_Amp_2IQi=Amp_all{2,ses}(Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))>=q1 & Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))<q2);
        FG_Amp_2IQi=Amp_all{3,ses}(Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))>=q1 & Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))<q2);
        
        [MI_SG_2IQi,MeanAmps_SG_2IQi]=ModIndex(Phase_2IQi, SG_Amp_2IQi, Ph_pos);
        [MI_FG_2IQi,MeanAmps_FG_2IQi]=ModIndex(Phase_2IQi, FG_Amp_2IQi, Ph_pos);
        
        Phase_3IQi=Phase_rad_all{1,ses}(Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))>=q2 & Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))<q3);
        SG_Amp_3IQi=Amp_all{2,ses}(Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))>=q2 & Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))<q3);
        FG_Amp_3IQi=Amp_all{3,ses}(Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))>=q2 & Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))<q3);
        
        [MI_SG_3IQi,MeanAmps_SG_3IQi]=ModIndex(Phase_3IQi, SG_Amp_3IQi, Ph_pos);
        [MI_FG_3IQi,MeanAmps_FG_3IQi]=ModIndex(Phase_3IQi, FG_Amp_3IQi, Ph_pos);

        Phase_4IQi=Phase_rad_all{1,ses}(Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))>=q3);
        SG_Amp_4IQi=Amp_all{2,ses}(Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))>=q3);
        FG_Amp_4IQi=Amp_all{3,ses}(Instpower_all{1,ses}(min_borders(minut)+1:min_borders(minut+1))>=q3);
        
        [MI_SG_4IQi,MeanAmps_SG_4IQi]=ModIndex(Phase_4IQi, SG_Amp_4IQi, Ph_pos);
        [MI_FG_4IQi,MeanAmps_FG_4IQi]=ModIndex(Phase_4IQi, FG_Amp_4IQi, Ph_pos);
        
        %now I collect all the outputs. First the mean gamma amplitudes for
        %each theta bin, for all 4 theta power IQis 
        Gamma_MeanAmp{1,1}(countmins,:)=MeanAmps_SG_1IQi;
        Gamma_MeanAmp{2,1}(countmins,:)=MeanAmps_FG_1IQi;
        Gamma_MeanAmp{1,2}(countmins,:)=MeanAmps_SG_2IQi;
        Gamma_MeanAmp{2,2}(countmins,:)=MeanAmps_FG_2IQi;
        Gamma_MeanAmp{1,3}(countmins,:)=MeanAmps_SG_3IQi;
        Gamma_MeanAmp{2,3}(countmins,:)=MeanAmps_FG_3IQi;
        Gamma_MeanAmp{1,4}(countmins,:)=MeanAmps_SG_4IQi;
        Gamma_MeanAmp{2,4}(countmins,:)=MeanAmps_FG_4IQi;
        
        SG_ModIndex(1,countmins)= MI_SG_1IQi;  
        FG_ModIndex(1,countmins)= MI_FG_1IQi;
        SG_ModIndex(2,countmins)= MI_SG_2IQi;
        FG_ModIndex(2,countmins)= MI_FG_2IQi;
        SG_ModIndex(3,countmins)= MI_SG_3IQi;
        FG_ModIndex(3,countmins)= MI_FG_3IQi;
        SG_ModIndex(4,countmins)= MI_SG_4IQi;
        FG_ModIndex(4,countmins)= MI_FG_4IQi;
      
  
        countmins=countmins+1;
          
    end 
end




%% Checking recurrence in PAC on a session


% Define the amplitude- and phase-frequencies
PhaseFreqVector=4:1:12; %this is the range of frequencies I'll be checking if their phase modulates the amplitude of others
AmpMGVector=60:2:100;
AmpSGVector=20:2:50;

PhaseFreq_BandWidth=2; %this the freqeuncy window size that'll be used in every iteration of the loop below
AmpFreq_BandWidth=10;


% define analyzed frequency ranges
PhaseFreqVectorPlot=2:2:50;
AmpFreqVectorPlot=10:5:200;
PhaseFreq_BandWidthPlot=4;
AmpFreq_BandWidthPlot=20;



%these next few loops (until line 757) are only to create comodulogram I
%can plot not for analysis. See after 756 for explanation
AmpTransformedPlot = cell(6,10); AmpTransformedPlot(:)={NaN(Fs*60,length(AmpFreqVectorPlot))}; 
PhaseFreqTransformedPlot = cell(6,10); PhaseFreqTransformedPlot(:)={NaN(Fs*60,length(PhaseFreqVectorPlot))};

for ses=1:sessions %for every session. That will also index the row in my cell array AmpFreqTransformed
    
    for minut=1:length(min_borders)-1 %this defines separate minutes within a whole session of Fs*600s = 2880000 samples    
        for ii=1:length(AmpFreqVectorPlot) 
  
            Af1 = AmpFreqVectorPlot(ii);
            Af2=Af1+AmpFreq_BandWidthPlot;
            AmpFreq=eegfilt((EGF_all_uV{ses}(min_borders(minut)+1:min_borders(minut+1)))',Fs,Af1,Af2); %I'm filtering here with this function and not the fftbandpass function because it only need and upper and lower threshold and it's fast. eegfilt needs row data input therefor I'm transposing
            AmpTransformedPlot{ses, minut}(:,ii) = abs(hilbert(AmpFreq')); % getting the amplitude envelope. I'm transposing AmpFreq back bcz I want it as a column and eegfilt gives it as a row 

        end 
    end
end

for ses=1:sessions 
    for minut=1:length(min_borders)-1
        for jj=1:length(PhaseFreqVectorPlot)
            
            Pf1 = PhaseFreqVectorPlot(jj);
            Pf2 = Pf1 + PhaseFreq_BandWidthPlot;
            PhaseFreq=eegfilt((EGF_all_uV{ses}(min_borders(minut)+1:min_borders(minut+1)))',Fs,Pf1,Pf2); % filtering 
            PhaseFreqTransformedPlot{ses, minut}(:,jj) = angle(hilbert(PhaseFreq')); % getting the phase time series
            
        end      
    end
end


ComodulogramPlot_all = cell(1,60); ComodulogramPlot_all(:)={NaN(length(PhaseFreqVectorPlot),length(AmpFreqVectorPlot))};
ComodulogramPlot_Zscored_all = cell(1,60); ComodulogramPlot_Zscored_all(:)={NaN(length(PhaseFreqVectorPlot),length(AmpFreqVectorPlot))};

minutcount=1;
for ses=1:sessions   
    disp(['Calculating wide frequency Phase Amplitude comodulograms for all minutes of session ', num2str(ses)]);
    for minut=1:length(min_borders)-1

        for jj=1:length(PhaseFreqVectorPlot) %for each frequency bin from the modulating frequency range
            for ii=1:length(AmpFreqVectorPlot) %check if it modulates the amplitude of all the frequency bins

                [MI,~]=ModIndex(PhaseFreqTransformedPlot{ses,minut}(:,jj), AmpTransformedPlot{ses,minut}(:,ii), Ph_pos);
                
                ComodulogramPlot_all{minutcount}(jj,ii)=MI;
            
            end 
        end
        %z-scoring the MI matrix to account for changes in the MI due to
        %electrode position
        ComodulogramPlot_Zscored_all{minutcount} = (ComodulogramPlot_all{minutcount}-mean(ComodulogramPlot_all{minutcount}(:)))./ std(ComodulogramPlot_all{minutcount}(:));
        minutcount=minutcount+1;
    end
end

% this is to average the comodulograms across all ten minutes of each
% session in an element-wise way. To do that I stack all 10 comodulograms
% of each session in a 3D matrix and then I average across the 3rd
% dimension.

ComodulogramPlot_3D=cell(1,6);
ComodulogramPlot_ses=cell(1,6);
ses=1;
for sesminstart=1:10:60
    ComodulogramPlot_3D{ses}=ComodulogramPlot_all{sesminstart};
    for sesmins=sesminstart+1:sesminstart+9
        ComodulogramPlot_3D{ses}=cat(3,ComodulogramPlot_3D{ses},ComodulogramPlot_all{sesmins});
    end
    ComodulogramPlot_ses{ses}=mean(ComodulogramPlot_3D{ses},3);
    ses=ses+1;
end



%Here I do what I did above but separetely for SG ang MG ranges. I'm filterring and creating amplitude and phase vectors for every minute and every session  

%Collecting the output of the loops
%The number of columns [length(AmpFreqVector) and length(PhaseFreqVector)] are different but every one of those columns is the same length (1min worth of data points)
%I'll then pass through the ModIndex function every column of
%AmpFreqTransformed against every colum of PhaseFreqTransformed (from the
%equivalent cells/mins) to create the comodulogram for every minute of the
%experiment. 
AmpSGTransformed = cell(6,10); AmpMSTransformed(:)={NaN(Fs*60,length(AmpSGVector))}; 
AmpMGTransformed = cell(6,10); AmpMGTransformed(:)={NaN(Fs*60,length(AmpMGVector))}; 
PhaseFreqTransformed = cell(6,10); PhaseFreqTransformed(:)={NaN(Fs*60,length(PhaseFreqVector))};

%what I'm doing with these loops is that I'm repeatedly filtering that eeg
%data (on a minute by minute basis) for different windows of frequency.
%These will essentially form the different bin in the comodulogram matrix

for ses=1:sessions %for every session. That will also index the row in my cell array AmpFreqTransformed
    
    for minut=1:length(min_borders)-1 %this defines separate minutes within a whole session of Fs*600s = 2880000 samples    
        for ii=1:length(AmpSGVector) 
  
            Af1 = AmpSGVector(ii);
            Af2=Af1+AmpFreq_BandWidth;
            AmpFreq=eegfilt((EGF_all_uV{ses}(min_borders(minut)+1:min_borders(minut+1)))',Fs,Af1,Af2); %I'm filtering here with this function and not the fftbandpass function because it only need and upper and lower threshold and it's fast. eegfilt needs row data input therefor I'm transposing
            AmpSGTransformed{ses, minut}(:,ii) = abs(hilbert(AmpFreq')); % getting the amplitude envelope. I'm transposing AmpFreq back bcz I want it as a column and eegfilt gives it as a row 

        end 
    end
end



for ses=1:sessions %for every session. That will also index the row in my cell array AmpFreqTransformed
    
    for minut=1:length(min_borders)-1 %this defines separate minutes within a whole session of Fs*600s = 2880000 samples    
        for ii=1:length(AmpMGVector) 
  
            Af1 = AmpMGVector(ii);
            Af2=Af1+AmpFreq_BandWidth;
            AmpFreq=eegfilt((EGF_all_uV{ses}(min_borders(minut)+1:min_borders(minut+1)))',Fs,Af1,Af2); %I'm filtering here with this function and not the fftbandpass function because it only need and upper and lower threshold and it's fast. eegfilt needs row data input therefor I'm transposing
            AmpMGTransformed{ses, minut}(:,ii) = abs(hilbert(AmpFreq')); % getting the amplitude envelope. I'm transposing AmpFreq back bcz I want it as a column and eegfilt gives it as a row 

        end 
    end
end

%repeat for the frequencies that have a modulating phase (usually the lower ones)
for ses=1:sessions 
    for minut=1:length(min_borders)-1
        for jj=1:length(PhaseFreqVector)
            
            Pf1 = PhaseFreqVector(jj);
            Pf2 = Pf1 + PhaseFreq_BandWidth;
            PhaseFreq=eegfilt((EGF_all_uV{ses}(min_borders(minut)+1:min_borders(minut+1)))',Fs,Pf1,Pf2); % filtering 
            PhaseFreqTransformed{ses, minut}(:,jj) = angle(hilbert(PhaseFreq')); % getting the phase time series
            
        end      
    end
end



ComodulogramMG_all = cell(1,60); ComodulogramMG_all(:)={NaN(length(PhaseFreqVector),length(AmpMGVector))}; %this has a single row with all the minute comodulograms to make it easier when I index this for reccurence
ComodulogramMG_Zscored_all = cell(1,60); ComodulogramMG_Zscored_all(:)={NaN(length(PhaseFreqVector),length(AmpMGVector))}; %this one will be used for the pearsons correlations across the different minutes

ComodulogramSG_all = cell(1,60); ComodulogramSG_all(:)={NaN(length(PhaseFreqVector),length(AmpSGVector))};
ComodulogramSG_Zscored_all = cell(1,60); ComodulogramSG_Zscored_all(:)={NaN(length(PhaseFreqVector),length(AmpSGVector))};

minutcount=1;
for ses=1:sessions   
    disp(['Calculating the SG-Theta Phase Amplitude comodulograms for all minutes of session ', num2str(ses)]);
    for minut=1:length(min_borders)-1

        for jj=1:length(PhaseFreqVector) %for each frequency bin from the modulating frequency range
            for ii=1:length(AmpSGVector) %check if it modulates the amplitude of all the frequency bins

                [MI,~]=ModIndex(PhaseFreqTransformed{ses,minut}(:,jj), AmpSGTransformed{ses,minut}(:,ii), Ph_pos);
                
                ComodulogramSG_all{minutcount}(jj,ii)=MI;
            
            end 
        end
        %z-scoring the MI matrix to account for changes in the MI due to
        %electrode position
        ComodulogramSG_Zscored_all{minutcount} = (ComodulogramSG_all{minutcount}-mean(ComodulogramSG_all{minutcount}(:)))./ std(ComodulogramSG_all{minutcount}(:));
        minutcount=minutcount+1;
    end
end


ComodulogramSG_3D=cell(1,6);
ComodulogramSG_ses=cell(1,6);
ses=1;
for sesminstart=1:10:60
    ComodulogramSG_3D{ses}=ComodulogramSG_all{sesminstart};
    for sesmins=sesminstart+1:sesminstart+9
        ComodulogramSG_3D{ses}=cat(3,ComodulogramSG_3D{ses},ComodulogramSG_all{sesmins});
    end
    ComodulogramSG_ses{ses}=mean(ComodulogramSG_3D{ses},3);
    ses=ses+1;
end



minutcount=1;
for ses=1:sessions   
    disp(['Calculating the MG-Theta Phase Amplitude comodulograms for all minutes of session ', num2str(ses)]);
    for minut=1:length(min_borders)-1

        for jj=1:length(PhaseFreqVector) %for each frequency bin from the modulating frequency range
            for ii=1:length(AmpMGVector) %check if it modulates the amplitude of all the frequency bins

                [MI,~]=ModIndex(PhaseFreqTransformed{ses,minut}(:,jj), AmpMGTransformed{ses,minut}(:,ii), Ph_pos);
                
                ComodulogramMG_all{minutcount}(jj,ii)=MI;
            
            end 
        end
        %z-scoring the MI matrix to account for changes in the MI due to
        %electrode position
        ComodulogramMG_Zscored_all{minutcount} = (ComodulogramMG_all{minutcount}-mean(ComodulogramMG_all{minutcount}(:)))./ std(ComodulogramMG_all{minutcount}(:));
        minutcount=minutcount+1;
    end
end


ComodulogramMG_3D=cell(1,6);
ComodulogramMG_ses=cell(1,6);
ses=1;
for sesminstart=1:10:60
    ComodulogramMG_3D{ses}=ComodulogramMG_all{sesminstart};
    for sesmins=sesminstart+1:sesminstart+9
        ComodulogramMG_3D{ses}=cat(3,ComodulogramMG_3D{ses},ComodulogramMG_all{sesmins});
    end
    ComodulogramMG_ses{ses}=mean(ComodulogramMG_3D{ses},3);
    ses=ses+1;
end










%% Checking recurrence in PAC on a session by session basis for all Theta power IQis
disp('Computing recurrence of PAC patterns across the experiment')

PAC_GammaTheta_Plot_minbymin_cor= NaN(length(ComodulogramPlot_Zscored_all),length(ComodulogramPlot_Zscored_all));

for corX=1:length(ComodulogramPlot_Zscored_all)
     for corY=1:length(ComodulogramPlot_Zscored_all)
         %if they're the same (autocorrelation) continue to leave NaNs
         %there. The reason I do that is to avoid having infin when I do
         %the fisher Z-transformation
         if corX==corY
             continue
         end
        
         cor_r= corr(ComodulogramPlot_Zscored_all{corX}(:),ComodulogramPlot_Zscored_all{corY}(:));
         PAC_GammaTheta_Plot_minbymin_cor(corX,corY)= 0.5*(log(1+cor_r) - log(1-cor_r)); %Fisher Z-Transformation of the pearson's corel coef to make it normally distributed
        
     end
 end    
    

PAC_GammaTheta_Plot_minbymin_cor=flipud(PAC_GammaTheta_Plot_minbymin_cor);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PAC_SG_cor_minbymin_cor= NaN(length(ComodulogramSG_Zscored_all),length(ComodulogramSG_Zscored_all));
PAC_MG_cor_minbymin_cor= NaN(length(ComodulogramMG_Zscored_all),length(ComodulogramMG_Zscored_all)); %the contain the autocor matrices for all mins

for corX=1:length(ComodulogramSG_Zscored_all)
     for corY=1:length(ComodulogramSG_Zscored_all)
         %if they're the same (autocorrelation) continue to leave NaNs
         %there. The reason I do that is to avoid having infin when I do
         %the fisher Z-transformation
         if corX==corY
             continue
         end     
         cor_r= corr(ComodulogramSG_Zscored_all{corX}(:),ComodulogramSG_Zscored_all{corY}(:));
         PAC_SG_cor_minbymin_cor(corX,corY)= 0.5*(log(1+cor_r) - log(1-cor_r)); %Fisher Z-Transformation of the pearson's corel coef to make it normally distributed
        
     end
 end    
    

PAC_SG_cor_minbymin_cor=flipud(PAC_SG_cor_minbymin_cor);


%%%%%%%%%%%%%%%

 for corX=1:length(ComodulogramMG_Zscored_all)
     for corY=1:length(ComodulogramMG_Zscored_all)
         %if they're the same (autocorrelation) continue to leave NaNs
         %there. The reason I do that is to avoid having infin when I do
         %the fisher Z-transformation
         if corX==corY
             continue
         end
        
         cor_r= corr(ComodulogramMG_Zscored_all{corX}(:),ComodulogramMG_Zscored_all{corY}(:));
         PAC_MG_cor_minbymin_cor(corX,corY)= 0.5*(log(1+cor_r) - log(1-cor_r)); %Fisher Z-Transformation of the pearson's corel coef to make it normally distributed
        
     end
 end    
    

PAC_MG_cor_minbymin_cor=flipud(PAC_MG_cor_minbymin_cor);





%% Section Four. Importing spike timings for the clusters indicated in the Clu2use xl file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

disp('Getting good cluster IDs from xl file...');
clu2use_cell=readcell('clu2use.xlsx'); %open the file which contains the usable clusters-needs to be created manually still
%clu2use_cell is a cell arrays and the first two columns have the rat ID
%and the date. I can use that in the exporting stage to make sure the the
%cell IDs are visible so I can import the data in the big matrix with the
%other information
Clu2useID=NaN(length(clu2use_cell),2);
for i=1:length(clu2use_cell)
    Clu2useID(i,1)=clu2use_cell{i,3}; %the 3rd column has the tetrode number 
    Clu2useID(i,2)=clu2use_cell{i,4}; %the fourth column has the cluster number for that tetrode
end

disp([num2str(length(Clu2useID)), ' good clusters in total in this session']);
Tetr= unique(Clu2useID(:,1)); %find the tetrodes that had good clusters
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
T4GoodClu=cell(1,length(Clu2useID)); %a cell array with the same dimensions as GoodCluindx_all
for ii=1:length(TtrID)
    GT=Feat_all{TtrID(ii)}(24,:); %getting the times for all the spikes in that tetrode
    T4GoodClu{ii}=GT(GoodCluindx_all{ii}); %getting only the spikes for the good clusters of that tetrode 
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

%% Section Five %%%% Bringing spikes and LFP together


%DOUBLECHECK that the inputs used here are the right ones!! AND ADD INFO
%FROM THE FIRST TWO COLUMNS OF clu2use_cell TO BE ABLE TO TRACK CLUSTERS

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
            SpkPhaseT_all{clu,ses}(spkT)= Phase_all{1,ses}(ix); %for theta
            
            SpkPhaseS_all{clu,ses}(spkT)= Phase_all{2,ses}(ix); %for slow gamma
            
            SpkPhaseF_all{clu,ses}(spkT)= Phase_all{3,ses}(ix); %for fast gamma
          
        end
    end
end


% Now I'll get the spikes that happend during gamma events. For theta is the entire session

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
        
        %this is to skip to the next iteration of the loop if a cluster has less than 12 spikes during gamma events (the number of phase bins) during a specific session
        %Here I skip the entire iteration (for all frequency bands) if the spikes are not enough for either slow or fast gamma...not very elegant
        if length(SPKPhaseS{clu,ses})<12 || length(SPKPhaseF{clu,ses})<12    
            continue
        end
        
        
        
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
             
        % Organise spike phases in a histogram and normalize spike count based on the number of spikes in the session/period of strong gamma. 
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
mkdir LFPAnalysisOutputnew %creates the folder in which the figures and xl files will be stored

baseFileName = 'PWRanalysis.xls'; %creating file that contains the data for the power septral density analysis. Section Two-line 317 onwards
fullFileName = fullfile('LFPAnalysisOutputnew', baseFileName);
csvwrite(fullfile('LFPAnalysisOutputnew','Vel_500ms.csv'),Velocity_cmpersec_500);
xlswrite(fullFileName ,vertcat(f_spect,Spectrogram_sessions_all(:,1:6)'),'Spctr_all');
xlswrite(fullFileName ,vertcat(f_spect,Spectrogram_sessions_all(:,7:12)'),'Spctr_Z_all');
xlswrite(fullFileName ,vertcat(f_spect,Spectrogram_log_sessions_all(:,1:6)'),'logSpctr_all');
xlswrite(fullFileName ,vertcat(f_spect,Spectrogram_log_sessions_all(:,7:12)'),'logSpctr_Z_all');
xlswrite(fullFileName ,Bands_all_ses,'Bands_all_ses'); 
xlswrite(fullFileName ,vertcat(f_spect,Spectrogram_move(:,1:6)'),'Spctr_mov_all');
xlswrite(fullFileName ,vertcat(f_spect,Spectrogram_move(:,7:12)'),'Spctr_mov_Z_all');
xlswrite(fullFileName ,vertcat(f_spect,Spectrogram_log_move(:,1:6)'),'logSpctr_mov_all');
xlswrite(fullFileName ,vertcat(f_spect,Spectrogram_log_move(:,7:12)'),'logSpctr_mov_Z_all');
xlswrite(fullFileName ,Bands_move_ses,'Bands_mov_ses'); 

%this is for individual bands theta, slow gamma and fast gamma very 2mins for uV EEG and z-scored. The first column is the average speed
%for that 2min period. This is without removing periods of imobility 
sesst=1;
for ses=1:sessions
    xlswrite(fullFileName ,Bands_2mins{ses},'PWR_bands_2min',['A', num2str(sesst)]);
    sesst=sesst+6;
end

indxl=2;
xlswrite(fullFileName,vbin_middle,'PWRsvsVel','A1')
xlswrite(fullFileName,vbin_middle,'PWRsvsVel','V1')
for bnd1=1:size(BandsPWR_Velbins,1) 
        
    xlswrite(fullFileName,BandsPWR_Velbins{bnd1,1},'PWRsvsVel',['A',num2str(indxl)])    
    xlswrite(fullFileName,BandsPWR_Velbins{bnd1,2},'PWRsvsVel',['V',num2str(indxl)])
      
    indxl=indxl+6;
    
end

Ph_bin_middle=10:20:350;
%This is for the oscillation phase asymmetries for every minute of the experiment_related to Kitanishi et al
%Figure 2H
xlswrite(fullFileName ,Ph_bin_middle,'ThetaPhasesymmetry','A1');
xlswrite(fullFileName ,Ph_bin_middle,'SgammaPhasesymmetry','A1');
xlswrite(fullFileName ,Ph_bin_middle,'MgammaPhasesymmetry','A1');
idxmin=2;
for ses=1:sessions
    
    xlswrite(fullFileName ,Phase_portion_all{1,ses},'ThetaPhasesymmetry',['A', num2str(idxmin)]);
    xlswrite(fullFileName ,Phase_portion_all{2,ses},'SgammaPhasesymmetry',['A', num2str(idxmin)]);
    xlswrite(fullFileName ,Phase_portion_all{3,ses},'MgammaPhasesymmetry',['A', num2str(idxmin)]);
    
    idxmin=idxmin+10;
end


%Now for the PAC (theta-gamma) related analysis.
baseFileName = 'PACrelatedanalysis.xls'; %creating file that contains the data for the Phase Amplitude coherence related analyses. Section Two-line 450 onwards
fullFileName = fullfile('LFPAnalysisOutputnew', baseFileName);

%this have the PAC MI means for all 4 theta IQis and all minutes of the experiment. 
%speed forthe individual theta power quantiles cannot be calculated because
%this is done using the instanteneous power (hilbert transform) and
%therefore I cannot reference the to the velocity which is calculated for
%every 500ms bin
xlswrite(fullFileName ,SG_ModIndex,'PAC_MI_min_SG');
xlswrite(fullFileName ,FG_ModIndex,'PAC_MI_min_MG');

%This is again for the modulation of gamma amplitude by theta phase. It's
%related to Tort et al 2009 PNAS Figure 4A
xlswrite(fullFileName ,Ph_bin_middle,'SGAmps','A1');
xlswrite(fullFileName ,Ph_bin_middle,'MGAmps','A1');
idxiqi=1;
for iqi=1:size(Gamma_MeanAmp,2)
    
    xlswrite(fullFileName ,Ph_bin_middle,'SGAmps',[xlscol(idxiqi),'1']);
    xlswrite(fullFileName ,Ph_bin_middle,'MGAmps',[xlscol(idxiqi),'1']);

    xlswrite(fullFileName ,Gamma_MeanAmp{1,iqi},'SGAmps',[xlscol(idxiqi),'2']);
    xlswrite(fullFileName ,Gamma_MeanAmp{2,iqi},'MGAmps',[xlscol(idxiqi),'2']);
    
    idxiqi=idxiqi+19;
end



%ADD THE ID OF THE CLUSTER HERE
%With this I collect the info I need so I can plot things like Kitanishi et al 2015 Neuron Figure 3B,D,E,F
Empty1=NaN(3,1);
for ses=1:size(ResulV_all,2)
    n=1;
    for clu=1:size(ResulV_all,1)
        
        xlswrite(fullFileName ,clu2use_cell{clu,3},['SpkPhaseLock_Session_', num2str(ses)],['A', num2str(n)]);
        xlswrite(fullFileName ,clu2use_cell{clu,4},['SpkPhaseLock_Session_', num2str(ses)],['B', num2str(n)]);
        
        if isempty(ResulV_all{clu,ses})
            
            xlswrite(fullFileName ,Empty1,['SpkPhaseLock_Session_', num2str(ses)],['C', num2str(n)]);
            xlswrite(fullFileName ,Empty1,['SpkPhaseLock_Session_', num2str(ses)],['D', num2str(n)]);
            xlswrite(fullFileName ,Empty1,['SpkPhaseLock_Session_', num2str(ses)],['E', num2str(n)]);
        else
            xlswrite(fullFileName ,ResulV_all{clu,ses},['SpkPhaseLock_Session_', num2str(ses)],['C', num2str(n)]);
            xlswrite(fullFileName ,MeanAng_all{clu,ses},['SpkPhaseLock_Session_', num2str(ses)],['D', num2str(n)]);
            xlswrite(fullFileName ,Raylp_all{clu,ses},['SpkPhaseLock_Session_', num2str(ses)],['E', num2str(n)]);
        end
   
        n=n+4;
    end
end


%With this I collect the info I need so I can plot things like Kitanishi et
%al 2015 Neuron Figure 3A. 
Empty2=NaN(3,12);
for ses=1:size(NormSpk_hist_all,2)
    xlswrite(fullFileName ,[x,x],['NormSpkCount_Session_', num2str(ses)],'C1');
    n=2;
    for clu=1:size(NormSpk_hist_all,1)
        
        xlswrite(fullFileName ,clu2use_cell{clu,3},['NormSpkCount_Session_', num2str(ses)],['A', num2str(n)]);
        xlswrite(fullFileName ,clu2use_cell{clu,4},['NormSpkCount_Session_', num2str(ses)],['B', num2str(n)]);
        
        if isempty(NormSpk_hist_all{clu,ses})
        
            xlswrite(fullFileName ,[Empty2,Empty2],['NormSpkCount_Session_', num2str(ses)],['C', num2str(n)]);
        else
            xlswrite(fullFileName ,[NormSpk_hist_all{clu,ses},NormSpk_hist_all{clu,ses}],['NormSpkCount_Session_', num2str(ses)],['C', num2str(n)]);
            
        end

        n=n+4;
    end
end

%Here I save the variabls that I didn't export in case Emma and Peter want
%something else in the analysis:

save([experiment,'extravariables','.mat'],'Bands_500bins','ComodulogramPlot_all','ComodulogramPlot_Zscored_all','ComodulogramPlot_ses','ComodulogramMG_all','ComodulogramMG_Zscored_all','ComodulogramSG_all','ComodulogramSG_Zscored_all','PAC_GammaTheta_Plot_minbymin_cor','PAC_SG_cor_minbymin_cor','PAC_MG_cor_minbymin_cor','-v7.3','-nocompression')  

%% Exporting figures

mkdir LFPAnalysisFiguresnew

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
    
    saveas(gcf,fullfile('LFPAnalysisFiguresnew',['Example traces for Session_',num2str(ses)]),'svg');
    close
end



for ses=1:sessions

    
    colormap jet
    contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpSGVector+AmpFreq_BandWidth/2,ComodulogramSG_ses{ses}',30,'lines','none');set(gcf,'Visible', 'off');
    set(gca,'fontsize',14)
    ylabel('Amplitude Frequency (Hz)')
    xlabel('Phase Frequency (Hz)')
    colorbar
    saveas(gcf,fullfile('LFPAnalysisFiguresnew',['PAC_SG-Th_Session ',num2str(ses)]),'svg');


   
    colormap jet
    contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpMGVector+AmpFreq_BandWidth/2,ComodulogramMG_ses{ses}',30,'lines','none');set(gcf,'Visible', 'off');
    set(gca,'fontsize',14)
    ylabel('Amplitude Frequency (Hz)')
    xlabel('Phase Frequency (Hz)')
    colorbar
    saveas(gcf,fullfile('LFPAnalysisFiguresnew',['PAC_MG-Th_Session ',num2str(ses)]),'svg');

end


colormap jet
imagesc(PAC_SG_cor_minbymin_cor);set(gcf,'Visible', 'off');
colorbar
saveas(gcf,fullfile('LFPAnalysisFiguresnew','PAC_SG_recurrence'),'svg');


colormap jet
imagesc(PAC_MG_cor_minbymin_cor);set(gcf,'Visible', 'off');
colorbar
saveas(gcf,fullfile('LFPAnalysisFiguresnew','PAC_MG_recurrence'),'svg');


colormap jet
imagesc(PAC_GammaTheta_Plot_minbymin_cor);set(gcf,'Visible', 'off');
colorbar
saveas(gcf,fullfile('LFPAnalysisFiguresnew','PAC_WideGamma_recurrence'),'svg');



disp('FINISHED!')
%toc;

%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%These are the some custom and modified functions that will be used in this analysis script.
%A few will be from the toolboxes FMAT, chronux, PAC and fieldtrip. I'll
%comment to indicate usage of toolbox functions

%% startup function for figures. This has been written by Fliss Inkpen, PhD

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

%% readEGF Provided by Takuma Kitanishi, PhD
% Custom function that reads high sampled eeg data and returns it in the array EEG together with the
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

%% EegMode Provided by Takuma Kitanishi, PhD

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
%    making the script chatty. Probably need to fix this in the future
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

%% fft function that is used in lines 449 and 510


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



%% eegfilt() -  (high|low|band)-iass filter data using two-way least-squares 
%              FIR filtering. Multiple data channels and epochs supported.
%              Requires the MATLAB Signal Processing Toolbox.
% Usage:
%  >> [smoothdata] = eegfilt(data,srate,locutoff,hicutoff);
%  >> [smoothdata,filtwts] = eegfilt(data,srate,locutoff,hicutoff, ...
%                                             epochframes,filtorder);
% Inputs:
%   data        = (channels,frames*epochs) data to filter
%   srate       = data sampling rate (Hz)
%   locutoff    = low-edge frequency in pass band (Hz)  {0 -> lowpass}
%   hicutoff    = high-edge frequency in pass band (Hz) {0 -> highpass}
%   epochframes = frames per epoch (filter each epoch separately {def/0: data is 1 epoch}
%   filtorder   = length of the filter in points {default 3*fix(srate/locutoff)}
%   revfilt     = [0|1] reverse filter (i.e. bandpass filter to notch filter). {0}
%
% Outputs:
%    smoothdata = smoothed data
%    filtwts    = filter coefficients [smoothdata <- filtfilt(filtwts,1,data)]
%
% See also: firls(), filtfilt()
%

% Author: Scott Makeig, Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 1997 

% Copyright (C) 4-22-97 from bandpass.m Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: eegfilt.m,v $
% Revision 1.19  2004/11/19 23:28:11  arno
% same
%
% Revision 1.18  2004/11/19 23:27:27  arno
% better error messsage
%
% Revision 1.17  2004/08/31 02:12:56  arno
% typo
%
% Revision 1.16  2004/02/12 23:09:19  scott
% same
%
% Revision 1.15  2004/02/12 23:08:02  scott
% text output edit
%
% Revision 1.14  2004/02/12 22:51:30  scott
% text output edits
%
% Revision 1.13  2003/07/28 17:38:53  arno
% updating error messages
%
% Revision 1.12  2003/07/20 19:17:07  scott
% added channels-processed info if epochs==1 (continuous data)
%
% Revision 1.11  2003/04/11 15:03:46  arno
% nothing
%
% Revision 1.10  2003/03/16 01:00:48  scott
% header and error msgs -sm
%
% Revision 1.9  2003/01/24 03:59:19  scott
% header edit -sm
%
% Revision 1.8  2003/01/24 00:23:33  arno
% print information about transition bands
%
% Revision 1.7  2003/01/23 23:53:25  arno
% change text
%
% Revision 1.6  2003/01/23 23:47:26  arno
% same
%
% Revision 1.5  2003/01/23 23:40:59  arno
% implementing notch filter
%
% Revision 1.4  2002/08/09 14:55:36  arno
% update transition band
%
% Revision 1.3  2002/08/09 02:06:27  arno
% updating transition
%
% Revision 1.2  2002/08/09 02:04:34  arno
% debugging
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 5-08-97 fixed frequency bound computation -sm
% 10-22-97 added MINFREQ tests -sm
% 12-05-00 added error() calls -sm
% 01-25-02 reformated help & license, added links -ad 

function [smoothdata,filtwts] = eegfilt(data,srate,locutoff,hicutoff,epochframes,filtorder, revfilt)

if nargin<4
    fprintf('');
    help eegfilt
    return
end

if ~exist('firls')
   error('*** eegfilt() requires the signal processing toolbox. ***');
end

[chans, frames] = size(data);
if chans > 1 && frames == 1
    help eegfilt
    error('input data should be a row vector.');
end
nyq            = srate*0.5;  % Nyquist frequency
%MINFREQ = 0.1/nyq;
MINFREQ = 0;

minfac         = 3;    % this many (lo)cutoff-freq cycles in filter 
min_filtorder  = 15;   % minimum filter length
trans          = 0.15; % fractional width of transition zones

if locutoff>0 && hicutoff > 0 && locutoff > hicutoff
    error('locutoff > hicutoff ???\n');
end
if locutoff < 0 || hicutoff < 0
   error('locutoff | hicutoff < 0 ???\n');
end

if locutoff>nyq
    error('Low cutoff frequency cannot be > srate/2');
end

if hicutoff>nyq
    error('High cutoff frequency cannot be > srate/2');
end

if nargin<6
   filtorder = 0;
end
if nargin<7
   revfilt = 0;
end

if isempty(filtorder) || filtorder==0
   if locutoff>0
     filtorder = minfac*fix(srate/locutoff);
   elseif hicutoff>0
     filtorder = minfac*fix(srate/hicutoff);
   end
     
   if filtorder < min_filtorder
        filtorder = min_filtorder;
    end
end

if nargin<5
	epochframes = 0;
end
if epochframes ==0
    epochframes = frames;    % default
end
epochs = fix(frames/epochframes);
if epochs*epochframes ~= frames
    error('epochframes does not divide frames.\n');
end

if filtorder*3 > epochframes   % Matlab filtfilt() restriction
    fprintf('eegfilt(): filter order is %d. ',filtorder);
    error('epochframes must be at least 3 times the filtorder.');
end
if (1+trans)*hicutoff/nyq > 1
	error('high cutoff frequency too close to Nyquist frequency');
end

if locutoff > 0 && hicutoff > 0    % bandpass filter
    if revfilt
%          fprintf('eegfilt() - performing %d-point notch filtering.\n',filtorder);
    else 
%          fprintf('eegfilt() - performing %d-point bandpass filtering.\n',filtorder);
    end 
%     fprintf('            If a message, ''Matrix is close to singular or badly scaled,'' appears,\n');
%     fprintf('            then Matlab has failed to design a good filter. As a workaround, \n');
%     fprintf('            for band-pass filtering, first highpass the data, then lowpass it.\n');

    f=[MINFREQ (1-trans)*locutoff/nyq locutoff/nyq hicutoff/nyq (1+trans)*hicutoff/nyq 1]; 
%     fprintf('eegfilt() - low transition band width is %1.1g Hz; high trans. band width, %1.1g Hz.\n',(f(3)-f(2))*srate, (f(5)-f(4))*srate/2);
    m=[0       0                      1            1            0                      0]; 
elseif locutoff > 0                % highpass filter
 if locutoff/nyq < MINFREQ
    error(sprintf('eegfilt() - highpass cutoff freq must be > %g Hz\n\n',MINFREQ*nyq));
 end
%  fprintf('eegfilt() - performing %d-point highpass filtering.\n',filtorder);
 f=[MINFREQ locutoff*(1-trans)/nyq locutoff/nyq 1]; 
%  fprintf('eegfilt() - highpass transition band width is %1.1g Hz.\n',(f(3)-f(2))*srate/2);
 m=[   0             0                   1      1];
elseif hicutoff > 0                %  lowpass filter
 if hicutoff/nyq < MINFREQ
    error(sprintf('eegfilt() - lowpass cutoff freq must be > %g Hz',MINFREQ*nyq));
 end
%  fprintf('eegfilt() - performing %d-point lowpass filtering.\n',filtorder);
 f=[MINFREQ hicutoff/nyq hicutoff*(1+trans)/nyq 1]; 
%  fprintf('eegfilt() - lowpass transition band width is %1.1g Hz.\n',(f(3)-f(2))*srate/2);
 m=[     1           1              0                 0];
else
    error('You must provide a non-0 low or high cut-off frequency');
end
if revfilt
    m = ~m;
end

filtwts = firls(filtorder,f,m); % get FIR filter coefficients

smoothdata = zeros(chans,frames);
for e = 1:epochs                % filter each epoch, channel 
    for c=1:chans
      smoothdata(c,(e-1)*epochframes+1:e*epochframes) ...
        = filtfilt(filtwts,1,data(c,(e-1)*epochframes+1:e*epochframes));
      if epochs == 1 
       if rem(c,20) ~= 0
%          fprintf('.');
       else 
%          fprintf('%d',c);
       end
      end
    end

end
% fprintf('\n');



end




%% Extract oscilation phase from the band pass filtered EEG signal by Hilbert transform. Troughs are set to 0/360 deg. Output is in deg
% September 2019 I added another output in rad because I need it for PAC. 
function [Phase, Phase_rad, Amp] = Osciphase(Eeg)

% Hilbert transform
Z = hilbert(Eeg);

% Wave amplitude
Amp = abs(Z);

% EEG phase in rad
Phase_rad = angle(Z);

% Rad to Degree (-180 to +180)
Phase = Phase_rad / pi *180;

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

%this number is based on the number of gamma events I usually get for 30-50Hz gamma for a 10min session. For 60-90Hz is about double.
%For theta peaks this will be met all the time, but I don't screen my
%spikes based on theta power peaks
if length(locs)<900 
    Th = mean(Power) + std(Power);
end

%retry to find peaks
[pks,locs] = findpeaks(Power);

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

%% [MI,MeanAmp]=ModIndex_v2(Phase, Amp, position) This is slightly modified for performance (for more details see Tort et al PNAS 2008, 2009 and J Neurophysiol 2010)
%
% Phase-amplitude cross-frequency coupling measure:
%
% Inputs:
% Phase = phase time series
% Amp = amplitude time series
% position = phase bins (left boundary)
%
% Outputs:
% MI = modulation index 
% MeanAmp = amplitude distribution over phase bins (non-normalized)
 
function [MI,MeanAmp]=ModIndex(Phase, Amp, position)

nbin=length(position);  
winsize = 2*pi/nbin;
 
% now we compute the mean amplitude in each phase:
MeanAmp=zeros(1,nbin); 
for j=1:nbin   
    MeanAmp(j) = mean(Amp(Phase <  position(j)+winsize & Phase >=  position(j)));
end
 
% the center of each bin (for plotting purposes) is position+winsize/2
 
% quantifying the amount of amp modulation by means of a
% normalized entropy index (Tort et al PNAS 2008):

MI=(log(nbin)-(-sum((MeanAmp/sum(MeanAmp)).*log((MeanAmp/sum(MeanAmp))))))/log(nbin);

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

%% MT spectrogram from FMAT toolbox. fliss Edit to add the log transformation

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

% Copyright (C) 2004-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [spectrogram, logTransformed,t_spect,f_spect] = MTSpectrogram(lfp,varargin)

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
if isempty(which('mtspecgramc'))
	error('This function requires the <a href="http://www.chronux.org">chronux</a> toolbox by P. Mitra, which does not appear to be installed on this system.');
end

% Check number of parameters
if nargin < 1 || mod(length(varargin),2) ~= 0
  error('Incorrect number of parameters (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
end

% Check parameter sizes
if size(lfp,2) ~= 1 && size(lfp,2) ~= 2
	error('Parameter ''lfp'' is not a vector or a Nx2 matrix (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin)
	if ~ischar(varargin{i})
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).']);
	end
	switch(lower(varargin{i}))
		case 'frequency'
			frequency = varargin{i+1};
			if ~isdscalar(frequency,'>0')
				error('Incorrect value for property ''frequency'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'range'
			range = varargin{i+1};
			if ~isdvector(range,'#2')
				error('Incorrect value for property ''range'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'window'
			window = varargin{i+1};
			if ~isdscalar(window,'>0')
				error('Incorrect value for property ''window'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'overlap'
			overlap = varargin{i+1};
			if ~isdscalar(overlap,'>0')
				error('Incorrect value for property ''overlap'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'step'
			step = varargin{i+1};
			if ~isdscalar(step,'>0')
				error('Incorrect value for property ''step'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'tapers'
			tapers = varargin{i+1};
			if ~isivector(tapers,'#2','>0')
				error('Incorrect value for property ''tapers'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'pad'
			pad = varargin{i+1};
			if ~isdscalar(pad,'>-1')
				error('Incorrect value for property ''pad'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'show'
			show = varargin{i+1};
			if ~isstring2(show,'on','off')
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help FindRipples">FindRipples</a>'' for details).');
			end
		case 'cutoffs'
			cutoffs = varargin{i+1};
			if ~isdvector(cutoffs,'#2','>=0')
				error('Incorrect value for property ''cutoffs'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
		case 'parent'
			parent = varargin{i+1};
			if ~ishandle(parent)
				error('Incorrect value for property ''parent'' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).');
			end
        otherwise
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help MTSpectrogram">MTSpectrogram</a>'' for details).']);
	end
end

% Determine LFP frequency if it is not already provided
if isempty(frequency)
	if size(lfp,2) == 2
		frequency = round(1/median(diff(lfp(:,1))));
	else
		frequency = f;
	end
end

% Determine step/overlap
if isempty(step)
	if isempty(overlap)
		overlap = window/2;
	end
else
	if isempty(overlap)
		overlap = window-step;
	else
		if overlap ~= window-step
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
    
logTransformed = log(spectrogram);
   
  
t_spect=t;
f_spect=f;
end

%% XLSCOL Convert Excel column letters to numbers or vice versa.
%   B = XLSCOL(A) takes input A, and converts to corresponding output B.
%   The input may be a number, a string, an array or matrix, an Excel
%   range, a cell, or a combination of each within a cell, including nested
%   cells and arrays. The output maintains the shape of the input and
%   attempts to "flatten" the cell to remove nesting.  Numbers and symbols
%   within strings or Excel ranges are ignored.
%
%   Examples
%   --------
%       xlscol(256)   % returns 'IV'
%
%       xlscol('IV')  % returns 256
%
%       xlscol([405 892])  % returns {'OO' 'AHH'}
%
%       xlscol('A1:IV65536')  % returns [1 256]
%
%       xlscol({8838 2430; 253 'XFD'}) % returns {'MAX' 'COL'; 'IS' 16384}
%
%       xlscol(xlscol({8838 2430; 253 'XFD'})) % returns same as input
%
%       b = xlscol({'A10' {'IV' 'ALL34:XFC66'} {'!@#$%^&*()'} '@#$' ...
%         {[2 3]} [5 7] 11})
%       % returns {1 [1x3 double] 'B' 'C' 'E' 'G' 'K'}
%       %   with b{2} = [256 1000 16383]
%
%   Notes
%   -----
%       CELLFUN and ARRAYFUN allow the program to recursively handle
%       multiple inputs.  An interesting side effect is that mixed input,
%       nested cells, and matrix shapes can be processed.
%
%   See also XLSREAD, XLSWRITE.
%
%   Version 1.1 - Kevin Crosby
% DATE      VER  NAME          DESCRIPTION
% 07-30-10  1.0  K. Crosby     First Release
% 08-02-10  1.1  K. Crosby     Vectorized loop for numerics.
% Contact: Kevin.L.Crosby@gmail.com

function [b] = xlscol(a)

base = 26;
if iscell(a)
  b = cellfun(@xlscol, a, 'UniformOutput', false); % handles mixed case too
elseif ischar(a)
  if ~isempty(strfind(a, ':')) % i.e. if is a range
    b = cellfun(@xlscol, regexp(a, ':', 'split'));
  else % if isempty(strfind(a, ':')) % i.e. if not a range
    b = a(isletter(a));        % get rid of numbers and symbols
    if isempty(b)
      b = {[]};
    else % if ~isempty(a);
      b = double(upper(b)) - 64; % convert ASCII to number from 1 to 26
      n = length(b);             % number of characters
      b = b * base.^((n-1):-1:0)';
    end % if isempty(a)
  end % if ~isempty(strfind(a, ':')) % i.e. if is a range
elseif isnumeric(a) && numel(a) ~= 1
  b = arrayfun(@xlscol, a, 'UniformOutput', false);
else % if isnumeric(a) && numel(a) == 1
  n = ceil(log(a)/log(base));  % estimate number of digits
  d = cumsum(base.^(0:n+1));   % offset
  n = find(a >= d, 1, 'last'); % actual number of digits
  d = d(n:-1:1);               % reverse and shorten
  r = mod(floor((a-d)./base.^(n-1:-1:0)), base) + 1;  % modulus
  b = char(r+64);  % convert number to ASCII
end % if iscell(a)
% attempt to "flatten" cell, by removing nesting
if iscell(b) && (iscell([b{:}]) || isnumeric([b{:}]))
  b = [b{:}];
end % if iscell(b) && (iscell([b{:}]) || isnumeric([ba{:}]))

end


