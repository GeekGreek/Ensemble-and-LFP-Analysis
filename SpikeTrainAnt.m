%Written by Antonis Asiminas Research Fellow SIDB, CDBS 
% May 2018
% 
% This is based on Neymotin et al J Neurosci Methods 2017. 
% This uses data from the Kluster-Neuroscope suite. 
% For more info on the file format architecture see http://neurosuite.sourceforge.net/formats.html
% It also uses an xl file with the identities of good quality clusters
% (based on isolation distnce and L-ratio). These were calculate
% previously using custom Matlab routines. For more info see Grieves, Wood, Dudchenko 2016 eLife
% Raw data were recorded in an Axona system. 

close all; %closes open figures and files
%% (Inputs)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sessions = input('How many sessions are you trying to analyse?');

%% (Section 1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elecs = [1 2 3 4 5 6 7 8]; %I don't really need this vector but it makes my life easier for the loops 
tic;
disp('Starting analysis...');
disp('Reading cluster identities...');
CluID_all = cell(1,8); %a 1x8 cell array with the cluster indentities for all 8 tetrodes
for i = 1:length(elecs)
        current_e = elecs(i);
clufile = ['merge.clu.' num2str(current_e)];
fid_clu = fopen(clufile,'r');
disp(sprintf('\t...getting cluster identities for tetrode %.f',i));
nclus = sscanf(fgetl(fid_clu),'%d');
cluID_elec = fscanf(fid_clu,'%f',[1 Inf]);
CluID_all{i}=cluID_elec;
fclose(fid_clu); %close the file at the end of each cycle to make the process faster
end
%% (Section 2)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reading spike features...');
Feat_all = cell(1,8); %a 1x8 cell array with the spike features for all 8 tetrodes
for k = 1:length(elecs)
        current_e = elecs(k);
fetfile = ['merge.fet.' num2str(current_e)];
fid_fet = fopen(fetfile, 'r');
disp(sprintf('\t...getting features for tetrode %.f',k));
nfeatures = sscanf(fgetl(fid_fet),'%d');
feat_elec = fscanf(fid_fet,'%f',[nfeatures Inf]);  
Feat_all{k}=feat_elec;
fclose(fid_fet);
end
%% (Section 3)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('Choosing usable clusters from xl file...');
Clu2useID = xlsread('Clu2use.xlsx'); %open the file which contains the usable clusters-needs to be created manually still
Tetr= unique(Clu2useID(:,1)); %find the tetrodes that had good clusters
disp(sprintf('\t%.f good clusters in total in this session',length(Clu2useID)));
GoodCluID_all = cell(1,8); %a 1x8 cell array with the good cluster indentities for all 8 tetrodes
for n=1:length(Tetr) %finding the usable clusters for each tetrode and creates a list for each tetrode
I=find(Clu2useID(:,1)==Tetr(n));
goodclusID = Clu2useID(I(1):I(end),2);
GoodCluID_all{Tetr(n)}= goodclusID;
end
%% (Section 4)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% (Section 5)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TtrID = Clu2useID(:,1);
T4GoodClu=cell(1,length(TtrID)); %a cell array with the same dimensions as GoodCluindx_all
for ii=1:length(TtrID)
    GT=Feat_all{TtrID(ii)}(24,:);
    t4goodclu=GT(GoodCluindx_all{ii});
    T4GoodClu{ii}=t4goodclu;
end

%T4GoodClu contains the timestamps for all the spikes of the good quality clusters 
%% (Section 6)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid_seg = fopen('merge.tseg', 'r'); %use the tseg file to see the size of the segments and the number of them
segs = fscanf(fid_seg,'%f',[2 Inf]);
T_segs =(segs(2,1:end))'; 
T_segs = [T_segs(1:end);T_segs(end)+T_segs(2)]*100000; %adds the start of first segment and the end of the last. For some reason this my be different for some data because it doesn't read the first zero
fclose(fid_seg);

bin_size = [2500 4000 12500 25000 100000 500000]; % these bins can change based on what I want. I may add an input argument later
Tbins_all = cell(1, (length(T_segs)-1)*length(bin_size)); 
h=0;
for f=1:length(T_segs)-1
   for b = bin_size
       bins=T_segs(f):b:T_segs(f+1);
       h=h+1;
       Tbins_all{h}=bins;
   end
end
%% (Section 7)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tbins_all contains all the tbin vectors for all sessions. All bin vectors for session 1, all bin vectors for sessions 2...
%i.e. For a session with 6 sessions and 6 time bins the array has 36 elements
h=0;
NoSpksInBins= cell(1,(length(T4GoodClu))*length(Tbins_all));
spksinbins =[];
for g=1:length(T4GoodClu) %for each good cluster
    disp(sprintf('\t...creating activity vectors for cluster %.f',g));
    for k=1:length(Tbins_all) %for each time bin vector
        for v=1:length(Tbins_all{k})-1 %for each time bin
            spkinbin=sum(T4GoodClu{g}>=Tbins_all{k}(v) & T4GoodClu{g}<=Tbins_all{k}(v+1)); %check how many spikes were there
            spksinbins(v)=spkinbin; %creat a vector for that cluster, for that bin vector
        end
    h=h+1;
        NoSpksInBins{h}=spksinbins;
        spksinbins =[];
    end
end

% NoSpksInBins a cell vector that contains the activity vectors for all clusters, all sessions and all bin sizes. All bins for all sessions for cluster 1, All bins for all sessions for cluster 2...and so on..
%i.e. For a session with 6 sessions and 6 time bins and 16 good clusters this array has 576 elements 
%% (Section 8)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%creation of the cell array Corr_pairs that contains the indexes of the all the
%cluster activity vectors (in NoSpksInBins) that need to be correlated in order to get the
%population correlation vectors(PCorrV) for all different bins and all different
%sessions. 
for v=1:length(Tbins_all)
k=v:length(Tbins_all):length(NoSpksInBins);
pout = nchoosek(k, 2);
Corr_pairs{v}=pout;
end
disp('Calculating Population correlation vectors...');
PCorr_Vs= cell(1,length(Corr_pairs)); %all PCorr vectors in a cell. It's structure is PCorrV for all(length(bin_size)) bin sizes for session 1, PCorrV for all bin sizes for session 2.....  
for v=1:length(Corr_pairs)
    for k=1:length(Corr_pairs{1}) %all the matrices in the cells of Corr_pairs are the same length since depend on the the number of clusters n*(n-1)/2
        Cl_corr = corr(NoSpksInBins{Corr_pairs{v}(k,1)}',NoSpksInBins{Corr_pairs{v}(k,2)}','Type','Kendall'); %gives Kendall's tau for every cluster combination in each session for the different time bins
        pcorr_v(k,1)=Cl_corr;
        pcorr_v(isnan(pcorr_v))=0; %replace NaN with 0. NaNs are generated when activity vectors are completely empty (only zeros)
    end
    PCorr_Vs{v}=pcorr_v;
end

%making the structure of PCorr_Vs more managable. Different rows are
%different bin sizes. It'll make it easy to plot later on
h=1;
for v=1:length(bin_size):length(PCorr_Vs)
    for k=1:length(bin_size)
    PCorr_Vs_Sessions{k,h}=PCorr_Vs{v+(k-1)}; %PCorr_Vs will be used to calculate Population Coordination in section 9.
    end
h=h+1;
end

%creating matrices for the different time (25 40 125 250 1000 5000ms) bins to make it easy to plot the PCorr vectors later on.

Mat_Sessions_bin25=zeros(length(PCorr_Vs_Sessions{1,1}),size(PCorr_Vs_Sessions,2));
for v=1:size(PCorr_Vs_Sessions,2)
 Mat_Sessions_bin25(:,v)=PCorr_Vs_Sessions{1,v};
end
Mat_Sessions_bin25=sortrows(Mat_Sessions_bin25,-1); %It sorts the PCorr values based on first column in descending order for better visual comparison when I generate the figures. The syntax of sortrows is different in future versions of MatLab

Mat_Sessions_bin40=zeros(length(PCorr_Vs_Sessions{1,1}),size(PCorr_Vs_Sessions,2));
for v=1:size(PCorr_Vs_Sessions,2)
 Mat_Sessions_bin40(:,v)=PCorr_Vs_Sessions{2,v};
end
Mat_Sessions_bin40=sortrows(Mat_Sessions_bin40,-1);

Mat_Sessions_bin125=zeros(length(PCorr_Vs_Sessions{1,1}),size(PCorr_Vs_Sessions,2));
for v=1:size(PCorr_Vs_Sessions,2)
 Mat_Sessions_bin125(:,v)=PCorr_Vs_Sessions{3,v};
end
Mat_Sessions_bin125=sortrows(Mat_Sessions_bin125,-1);

Mat_Sessions_bin250=zeros(length(PCorr_Vs_Sessions{1,1}),size(PCorr_Vs_Sessions,2));
for v=1:size(PCorr_Vs_Sessions,2)
 Mat_Sessions_bin250(:,v)=PCorr_Vs_Sessions{4,v};
end
Mat_Sessions_bin250=sortrows(Mat_Sessions_bin250,-1);

Mat_Sessions_bin1s=zeros(length(PCorr_Vs_Sessions{1,1}),size(PCorr_Vs_Sessions,2));
for v=1:size(PCorr_Vs_Sessions,2)
 Mat_Sessions_bin1s(:,v)=PCorr_Vs_Sessions{5,v};
end
Mat_Sessions_bin1s=sortrows(Mat_Sessions_bin1s,-1);

Mat_Sessions_bin5s=zeros(length(PCorr_Vs_Sessions{1,1}),size(PCorr_Vs_Sessions,2));
for v=1:size(PCorr_Vs_Sessions,2)
 Mat_Sessions_bin5s(:,v)=PCorr_Vs_Sessions{6,v};
end
Mat_Sessions_bin5s=sortrows(Mat_Sessions_bin5s,-1);

%% (Section 9)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating Population Coordination...');
session_pairs=nchoosek(1:length(bin_size):length(PCorr_Vs),2); %creating a matrix containing all the different pairs of sessions that need to be correlated to get the coordination of population.
session_pairs_PCo=nchoosek(1:length(T_segs)-1,2); %creates a matrix with the first two columns showing which sessions are being correlated to be used in the loop
for v=1:length(bin_size)
    for k=1:length(session_pairs)
        [PCo,p]=corr(PCorr_Vs{session_pairs(k,1)+(v-1)},PCorr_Vs{session_pairs(k,2)+(v-1)},'Type','Pearson');
        session_pairs_PCo(k,3)=PCo;
        session_pairs_PCo(k,4)=p;  
    end
    PCo_Allbins{v}=session_pairs_PCo;
end

%PCo_Allbins contains as many cells as the length of bin_size. Each cell
%contains a matrix in which the first 2 column show which sessions are
%being correlated (Pearson), the third column has the value r and the
%fourth column has the p value. It maybe the case that for this step I
%should put together all the pairs for all the rats and then correlate
%between sessions.

%% (Section 10)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%In this section I'll revisit section 6 and instead of splitting the entire
%recording to the number of sessions (4 or 6 in my case) I'll split each
%session in half and calculate the population correlation vectors (section 8) and then the Population Coordination (section 9) between the two half of each session.
    
T_segs_half =zeros((length(T_segs)-1)*2+1,1); % creates a column vector with the start and end of all sessions*2 epochs (sessions*2 because we have split the sessions into half)
    for n=2:length(T_segs)
        Session_halfn = (T_segs(n)-T_segs(n-1))/2;
        T_segs_half((n-1)*2)= T_segs(n-1)+Session_halfn;
        T_segs_half((n-1)*2+1)= T_segs(n); 
    end

%This creates exactly the same as in section 6 but with double the epochs
%because we split the epochs in half
Tbins_all_half = cell(1, (length(T_segs)-1)*length(bin_size)); 
h=0;
for f=1:length(T_segs_half)-1
   for b = bin_size
       bins=[T_segs_half(f):b:T_segs_half(f+1)];
       h=h+1;
       Tbins_all_half{h}=bins;
   end
end

%This creates exactly the same as in section 7 but with double the epochs
%because we split the epochs in half
h=0;
NoSpksInBins_half= cell(1,(length(T4GoodClu))*length(Tbins_all_half));
spksinbins =[];
for g=1:length(T4GoodClu) %for each good cluster
    disp(sprintf('\t...creating activity vectors for cluster %.f',g));
    for k=1:length(Tbins_all_half) %for each time bin vector
        for v=1:length(Tbins_all_half{k})-1 %for each time bin
            spkinbin=sum(T4GoodClu{g}>=Tbins_all_half{k}(v) & T4GoodClu{g}<=Tbins_all_half{k}(v+1)); %check how many spikes were there
            spksinbins(v)=spkinbin; %creat a vector for that cluster, for that bin vector
        end
        h=h+1;
        NoSpksInBins_half{h}=spksinbins;
        spksinbins =[];
    end
end

for v=1:length(Tbins_all_half)
k=v:length(Tbins_all_half):length(NoSpksInBins_half);
pout = nchoosek(k, 2);
Corr_pairs_half{v}=pout;
end
disp('Calculating Population correlation vectors for halves of sessions...');
PCorr_Vs_half= cell(1,length(Corr_pairs_half)); %all PCorr vectors in a cell. It's structure is PCorrV for all(length(bin_size)) bin sizes for session 1 first half, PCorrV for all bin sizes for session 1 second half.....  
for v=1:length(Corr_pairs_half)
    for k=1:length(Corr_pairs_half{1}) %all the matrices in the cells of Corr_pairs are the same length since depend on the the number of clusters n*(n-1)/2
        Cl_corr = corr(NoSpksInBins_half{Corr_pairs_half{v}(k,1)}',NoSpksInBins_half{Corr_pairs_half{v}(k,2)}','Type','Kendall'); %gives Kendall's tau for every cluster combination in each session half for the different time bins
        pcorr_vh(k,1)=Cl_corr;
        pcorr_vh(isnan(pcorr_vh))=0; %replace NaN with 0. NaNs are generated when activity vectors are completely empty (only zeros)
    end
    PCorr_Vs_half{v}=pcorr_vh;
end

%making the structure of PCorr_Vs_half more managable. Different row are
%different bin sizes. It'll make it easy to plot later on
h=1;
for v=1:length(bin_size):length(PCorr_Vs_half)
    for k=1:length(bin_size)
    PCorr_Vs_half_Sessions{k,h}=PCorr_Vs_half{v+(k-1)}; %PCorr_Vs_half will be used to calculate Population Coordination in section 9.
    end
h=h+1;
end

%creating matrices for the different time (25 40 125 250 1000 5000ms) bins to make it easy to plot the PCorr vectors later on.

Mat_half_Sessions_bin25=zeros(length(PCorr_Vs_half_Sessions{1,1}),size(PCorr_Vs_half_Sessions,2));
for v=1:size(PCorr_Vs_half_Sessions,2)
 Mat_half_Sessions_bin25(:,v)=PCorr_Vs_half_Sessions{1,v};
end

Mat_half_Sessions_bin40=zeros(length(PCorr_Vs_half_Sessions{1,1}),size(PCorr_Vs_half_Sessions,2));
for v=1:size(PCorr_Vs_half_Sessions,2)
 Mat_half_Sessions_bin40(:,v)=PCorr_Vs_half_Sessions{2,v};
end

Mat_half_Sessions_bin125=zeros(length(PCorr_Vs_half_Sessions{1,1}),size(PCorr_Vs_half_Sessions,2));
for v=1:size(PCorr_Vs_half_Sessions,2)
 Mat_half_Sessions_bin125(:,v)=PCorr_Vs_half_Sessions{3,v};
end

Mat_half_Sessions_bin250=zeros(length(PCorr_Vs_half_Sessions{1,1}),size(PCorr_Vs_half_Sessions,2));
for v=1:size(PCorr_Vs_half_Sessions,2)
 Mat_half_Sessions_bin250(:,v)=PCorr_Vs_half_Sessions{4,v};
end

Mat_half_Sessions_bin1s=zeros(length(PCorr_Vs_half_Sessions{1,1}),size(PCorr_Vs_half_Sessions,2));
for v=1:size(PCorr_Vs_half_Sessions,2)
 Mat_half_Sessions_bin1s(:,v)=PCorr_Vs_half_Sessions{5,v};
end

Mat_half_Sessions_bin5s=zeros(length(PCorr_Vs_half_Sessions{1,1}),size(PCorr_Vs_half_Sessions,2));
for v=1:size(PCorr_Vs_half_Sessions,2)
 Mat_half_Sessions_bin5s(:,v)=PCorr_Vs_half_Sessions{6,v};
end


%Calculating the Population coordination between the two half of each
%session only
disp('Calculating Population Coordination...');
index_half_sessions=1:length(bin_size):length(Tbins_all_half); 
session_halfpairs(:,1)=index_half_sessions(1:2:end);
session_halfpairs(:,2)=index_half_sessions(2:2:end); %creating a matrix containing all the different pairs of half sessions that need to be correlated to get the coordination of population.
session_half_PCo(:,1) =1:length(T_segs)-1; %creates the headers (first column) for the correlations calculated in the loop
for v=1:length(bin_size) 
    for k=1:length(session_halfpairs)
        [PCo,p]=corr(PCorr_Vs_half{session_halfpairs(k,1)+(v-1)},PCorr_Vs_half{session_halfpairs(k,2)+(v-1)},'Type','Pearson');
        session_half_PCo(k,2)=PCo;
        session_half_PCo(k,3)=p;  
    end
    PCo_HalfAll{v}=session_half_PCo; %the cells are as many bin sizes we have and inside the cells are matrices with the r and p values for the population coordination between 1st and 2nd half of each session
end

%% (Section 11)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%In this section I repeat all steps in section 10 but now I create 1 min
%epochs to create a big matrix for Population coordination

Minute = T_segs(2)/10; %every session is presumably the same length so this gives the length of 1 min from the first session
T_segs_min =[0:Minute:T_segs(end)]'; % creates a column vector with the start and end of all sessions*10 minutes 

%This creates exactly the same as in section 6 but with 10 times the epochs
%because we divided recording in minutes
Tbins_all_min = cell(1, (length(T_segs)-1)*length(bin_size)); 
h=0;
for f=1:length(T_segs_min)-1
   for b = bin_size
       bins=[T_segs_min(f):b:T_segs_min(f+1)];
       h=h+1;
       Tbins_all_min{h}=bins;
   end
end

%This creates exactly the same as in section 6 but with 10 times the epochs
%because we divided recording in minutes
h=0;
NoSpksInBins_min= cell(1,(length(T4GoodClu))*length(Tbins_all_min));
spksinbins =[];
for g=1:length(T4GoodClu) %for each good cluster
    disp(sprintf('\t...creating activity vectors for cluster %.f',g));
    for k=1:length(Tbins_all_min) %for each time bin vector
        for v=1:length(Tbins_all_min{k})-1 %for each time bin
            spkinbin=sum(T4GoodClu{g}>=Tbins_all_min{k}(v) & T4GoodClu{g}<=Tbins_all_min{k}(v+1)); %check how many spikes were there
            spksinbins(v)=spkinbin; %creat a vector for that cluster, for that bin vector
        end
        h=h+1;
        NoSpksInBins_min{h}=spksinbins;
        spksinbins =[];
    end
end

for v=1:length(Tbins_all_min)
k=v:length(Tbins_all_min):length(NoSpksInBins_min);
pout = nchoosek(k, 2);
Corr_pairs_min{v}=pout;
end
disp('Calculating Population correlation vectors for every minute of the sessions...');
PCorr_Vs_min= cell(1,length(Corr_pairs_min)); %all PCorr vectors in a cell. It's structure is PCorrV for all(length(bin_size)) bin sizes for 1st min of recording, PCorrV for all bin sizes for 2nd min of recording.....  
pcorr_vm=zeros(length(TtrID)*(length(TtrID)-1)/2, 1);
for v=1:length(Corr_pairs_min)
    for k=1:length(Corr_pairs_min{1}) %all the matrices in the cells of Corr_pairs are the same length since depend on the the number of clusters n*(n-1)/2
        Cl_corr = corr(NoSpksInBins_min{Corr_pairs_min{v}(k,1)}',NoSpksInBins_min{Corr_pairs_min{v}(k,2)}','Type','Kendall'); %gives Kendall's tau for every cluster combination in each min of recording for the different time bins
        pcorr_vm(k,1)=Cl_corr;
        pcorr_vm(isnan(pcorr_vm))=0; %replace NaN with 0. NaNs are generated when activity vectors are completely empty (only zeros)
    end
    PCorr_Vs_min{v}=pcorr_vm;
end

%making the structure of PCorr_Vs_min more managable. Different row are
%different bin sizes
h=1;
PCorr_Vs_Min=cell(sessions,length(T_segs_min)-1); 
for v=1:length(bin_size):length(PCorr_Vs_min)
    for k=1:length(bin_size)
    PCorr_Vs_Min{k,h}=PCorr_Vs_min{v+(k-1)};
    end
h=h+1;
end

%Calculating the Population coordination across the recording
disp('Calculating Population Coordination...');
%creating matrices for the different time (25 40 125 250 1000 5000ms) bins so I can run
%cross-correlation analysis for each of them and make it easy to plot the
%matrices later on.
Mat_bin25=zeros(length(PCorr_Vs_Min{1,1}),length(PCorr_Vs_Min));
for v=1:size(PCorr_Vs_Min,2)
 Mat_bin25(:,v)=PCorr_Vs_Min{1,v};
end
%Creating PCo matrix for 25ms bin
disp('Creating PCo matrix for 25ms bin...');
CorrMat25=zeros(length(T_segs_min)-1);
for m1 = 1:size(Mat_bin25,2) % Create correlations for each epoch
 for m2 = 1:size(Mat_bin25,2) % Correlate against each epoch
  CorrMat25(m2,m1) = corr(Mat_bin25(:,m2),Mat_bin25(:,m1));
 end
end
CorrMat25=flipud(CorrMat25);

Mat_bin40=zeros(length(PCorr_Vs_Min{1,1}),length(PCorr_Vs_Min));
for v=1:size(PCorr_Vs_Min,2)
 Mat_bin40(:,v)=PCorr_Vs_Min{2,v};
end
%creating PCo matrix for 40ms bin
disp('Creating PCo matrix for 40ms bin...');
CorrMat40=zeros(length(T_segs_min)-1);
for m1 = 1:size(Mat_bin40,2) % Create correlations for each epoch
 for m2 = 1:size(Mat_bin40,2) % Correlate against each epoch
  CorrMat40(m2,m1) = corr(Mat_bin40(:,m2),Mat_bin40(:,m1));
 end
end
CorrMat40=flipud(CorrMat40);

Mat_bin125=zeros(length(PCorr_Vs_Min{1,1}),length(PCorr_Vs_Min));
for v=1:size(PCorr_Vs_Min,2)
 Mat_bin125(:,v)=PCorr_Vs_Min{3,v};
end
%creating PCo matrix for 125ms bin
disp('Creating PCo matrix for 125ms bin...');
CorrMat125=zeros(length(T_segs_min)-1);
for m1 = 1:size(Mat_bin125,2) % Create correlations for each epoch
 for m2 = 1:size(Mat_bin125,2) % Correlate against each epoch
  CorrMat125(m2,m1) = corr(Mat_bin125(:,m2),Mat_bin125(:,m1));
 end
end
CorrMat125=flipud(CorrMat125);

Mat_bin250=zeros(length(PCorr_Vs_Min{1,1}),length(PCorr_Vs_Min));
for v=1:size(PCorr_Vs_Min,2)
 Mat_bin250(:,v)=PCorr_Vs_Min{4,v};
end
%creating PCo matrix for 250ms bin
disp('Creating PCo matrix for 250ms bin...');
CorrMat250=zeros(length(T_segs_min)-1);
for m1 = 1:size(Mat_bin250,2) % Create correlations for each epoch
 for m2 = 1:size(Mat_bin250,2) % Correlate against each epoch
  CorrMat250(m2,m1) = corr(Mat_bin250(:,m2),Mat_bin250(:,m1));
 end
end
CorrMat250=flipud(CorrMat250);

Mat_bin1s=zeros(length(PCorr_Vs_Min{1,1}),length(PCorr_Vs_Min));
for v=1:size(PCorr_Vs_Min,2)
 Mat_bin1s(:,v)=PCorr_Vs_Min{5,v};
end
%creating PCo matrix for 1s bin
disp('Creating PCo matrix for 1s bin...');
CorrMat1s=zeros(length(T_segs_min)-1);
for m1 = 1:size(Mat_bin1s,2) % Create correlations for each epoch
 for m2 = 1:size(Mat_bin1s,2) % Correlate against each epoch
  CorrMat1s(m2,m1) = corr(Mat_bin1s(:,m2),Mat_bin1s(:,m1));
 end
end
CorrMat1s=flipud(CorrMat1s);

Mat_bin5s=zeros(length(PCorr_Vs_Min{1,1}),length(PCorr_Vs_Min));
for v=1:size(PCorr_Vs_Min,2)
 Mat_bin5s(:,v)=PCorr_Vs_Min{6,v};
end
%creating PCo matrix for 5s bin
disp('Creating PCo matrix for 5s bin...');
CorrMat5s=zeros(length(T_segs_min)-1);
for m1 = 1:size(Mat_bin5s,2) % Create correlations for each epoch
 for m2 = 1:size(Mat_bin5s,2) % Correlate against each epoch
  CorrMat5s(m2,m1) = corr(Mat_bin5s(:,m2),Mat_bin5s(:,m1));
 end
end
CorrMat5s=flipud(CorrMat5s);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creating output files

mkdir SpikeTrainOutput %creates the folder in which the figures and xl files will be stored

baseFileName = 'Entire sessions.xls'; %creating file that contains the data for entire sessions and PCo between all sessions 
fullFileName = fullfile('SpikeTrainOutput', baseFileName);
xlswrite(fullFileName ,Mat_Sessions_bin25,'PCorr25ms');
xlswrite(fullFileName ,PCo_Allbins{1},'PCo25ms');

xlswrite(fullFileName ,Mat_Sessions_bin40,'PCorr40ms');
xlswrite(fullFileName ,PCo_Allbins{2},'PCo40ms');

xlswrite(fullFileName ,Mat_Sessions_bin125,'PCorr125ms');
xlswrite(fullFileName ,PCo_Allbins{3},'PCo125ms');

xlswrite(fullFileName ,Mat_Sessions_bin250,'PCorr250ms');
xlswrite(fullFileName ,PCo_Allbins{4},'PCo250ms');

xlswrite(fullFileName ,Mat_Sessions_bin1s,'PCorr1s');
xlswrite(fullFileName ,PCo_Allbins{5},'PCo1s');

xlswrite(fullFileName ,Mat_Sessions_bin5s,'PCorr5s');
xlswrite(fullFileName ,PCo_Allbins{6},'PCo5s');


baseFileName = 'Half sessions.xls'; %creating file that contains the data for half sessions and PCo between first and second half for each session 
fullFileName = fullfile('SpikeTrainOutput', baseFileName);
xlswrite(fullFileName ,Mat_half_Sessions_bin25,'PCorr25ms');
xlswrite(fullFileName ,PCo_HalfAll{1},'PCo25ms');

xlswrite(fullFileName ,Mat_half_Sessions_bin40,'PCorr40ms');
xlswrite(fullFileName ,PCo_HalfAll{2},'PCo40ms');

xlswrite(fullFileName ,Mat_half_Sessions_bin125,'PCorr125ms');
xlswrite(fullFileName ,PCo_HalfAll{3},'PCo125ms');

xlswrite(fullFileName ,Mat_half_Sessions_bin250,'PCorr250ms');
xlswrite(fullFileName ,PCo_HalfAll{4},'PCo250ms');

xlswrite(fullFileName ,Mat_half_Sessions_bin1s,'PCorr1s');
xlswrite(fullFileName ,PCo_HalfAll{5},'PCo1s');

xlswrite(fullFileName ,Mat_half_Sessions_bin5s,'PCorr5s');
xlswrite(fullFileName ,PCo_HalfAll{6},'PCo5s');

%creating file that contains the data for minbimin PCos_In case I want to plot them using a different software or a different script later_I'll need the values to get max and min to add in my figures' colourbars
baseFileName = 'MinbyMin.xls'; 
fullFileName = fullfile('SpikeTrainOutput', baseFileName);
xlswrite(fullFileName ,CorrMat25,'PCo25ms');

xlswrite(fullFileName ,CorrMat40,'PCo40ms');

xlswrite(fullFileName ,CorrMat125,'PCo125ms');

xlswrite(fullFileName ,CorrMat250,'PCo250ms');

xlswrite(fullFileName ,CorrMat1s,'PCo1s');

xlswrite(fullFileName ,CorrMat5s,'PCo5s');

%% Plotting PCorr vectors for entire sessions
format = 'png';
mkdir FiguresSpikes
%for 25ms
for k=1:sessions
    h=imagesc(Mat_Sessions_bin25(:,k));set(gcf,'Visible', 'off'); 
    colorbar 
    caxis([min(Mat_Sessions_bin25(:)) max(Mat_Sessions_bin25(:))]);
    pngFileName = sprintf('Entire_session_%dbin_25ms', k);
    saveas(gcf,fullfile('FiguresSpikes',pngFileName),format);
end

%for 40ms
for k=1:sessions
    h=imagesc(Mat_Sessions_bin40(:,k));set(gcf,'Visible', 'off'); 
    colorbar 
    caxis([min(Mat_Sessions_bin40(:)) max(Mat_Sessions_bin40(:))]);
    pngFileName = sprintf('Entire_session_%dbin_40ms', k);
    saveas(gcf,fullfile('FiguresSpikes',pngFileName),format);
end

%for 125ms
for k=1:sessions
    h=imagesc(Mat_Sessions_bin125(:,k));set(gcf,'Visible', 'off'); 
    colorbar 
    caxis([min(Mat_Sessions_bin125(:)) max(Mat_Sessions_bin125(:))]);
    pngFileName = sprintf('Entire_session_%dbin_125ms', k);
    saveas(gcf,fullfile('FiguresSpikes',pngFileName),format);
end

%for 250ms
for k=1:sessions
    h=imagesc(Mat_Sessions_bin250(:,k));set(gcf,'Visible', 'off'); 
    colorbar 
    caxis([min(Mat_Sessions_bin250(:)) max(Mat_Sessions_bin250(:))]);
    pngFileName = sprintf('Entire_session_%dbin_250ms', k);
    saveas(gcf,fullfile('FiguresSpikes',pngFileName),format);
end

%for 1s
for k=1:sessions
    h=imagesc(Mat_Sessions_bin1s(:,k));set(gcf,'Visible', 'off'); 
    colorbar 
    caxis([min(Mat_Sessions_bin1s(:)) max(Mat_Sessions_bin1s(:))]);
    pngFileName = sprintf('Entire_session_%dbin_2s', k);
    saveas(gcf,fullfile('FiguresSpikes',pngFileName),format);
end

%for 5s
for k=1:sessions
    h=imagesc(Mat_Sessions_bin5s(:,k));set(gcf,'Visible', 'off'); 
    colorbar 
    caxis([min(Mat_Sessions_bin5s(:)) max(Mat_Sessions_bin5s(:))]);
    pngFileName = sprintf('Entire_session_%dbin_5s', k);
    saveas(gcf,fullfile('FiguresSpikes',pngFileName),format);
end

%% Plotting PCorr vectors for hafves of sessions

%for 25ms
for k=1:2:length(T_segs_half)-1
    TempMat=Mat_half_Sessions_bin25(:,k:k+1);
    TempMat=sortrows(TempMat,-1);
    for n=1:size(TempMat,2)
        h=imagesc(TempMat(:,n));set(gcf,'Visible', 'off');  
        colorbar
        caxis([min(TempMat(:)) max(TempMat(:))]);
        pngFileName = sprintf('Session_%d_%dbin_25ms', (k+1)/2, n);
        saveas(gcf,fullfile('FiguresSpikes',pngFileName),format);
    end
end

%for 40ms
for k=1:2:length(T_segs_half)-1
    TempMat=Mat_half_Sessions_bin40(:,k:k+1);
    TempMat=sortrows(TempMat,-1);
    for n=1:size(TempMat,2)
        h=imagesc(TempMat(:,n));set(gcf,'Visible', 'off');  
        colorbar
        caxis([min(TempMat(:)) max(TempMat(:))]);
        pngFileName = sprintf('Session_%d_%dbin_40ms', (k+1)/2, n);
        saveas(gcf,fullfile('FiguresSpikes',pngFileName),format);
    end
end

%for 125ms
for k=1:2:length(T_segs_half)-1
    TempMat=Mat_half_Sessions_bin125(:,k:k+1);
    TempMat=sortrows(TempMat,-1);
    for n=1:size(TempMat,2)
        h=imagesc(TempMat(:,n));set(gcf,'Visible', 'off');  
        colorbar
        caxis([min(TempMat(:)) max(TempMat(:))]);
        pngFileName = sprintf('Session_%d_%dbin_125ms', (k+1)/2, n);
        saveas(gcf,fullfile('FiguresSpikes',pngFileName),format);
    end
end

%for 250ms
for k=1:2:length(T_segs_half)-1
    TempMat=Mat_half_Sessions_bin250(:,k:k+1);
    TempMat=sortrows(TempMat,-1);
    for n=1:size(TempMat,2)
        h=imagesc(TempMat(:,n));set(gcf,'Visible', 'off');  
        colorbar
        caxis([min(TempMat(:)) max(TempMat(:))]);
        pngFileName = sprintf('Session_%d_%dbin_250ms', (k+1)/2, n);
        saveas(gcf,fullfile('FiguresSpikes',pngFileName),format);
    end
end

%for 1s
for k=1:2:length(T_segs_half)-1
    TempMat=Mat_half_Sessions_bin1s(:,k:k+1);
    TempMat=sortrows(TempMat,-1);
    for n=1:size(TempMat,2)
        h=imagesc(TempMat(:,n));set(gcf,'Visible', 'off');  
        colorbar
        caxis([min(TempMat(:)) max(TempMat(:))]);
        pngFileName = sprintf('Session_%d_%dbin_1s', (k+1)/2, n);
        saveas(gcf,fullfile('FiguresSpikes',pngFileName),format);
    end
end

%for 5s
for k=1:2:length(T_segs_half)-1
    TempMat=Mat_half_Sessions_bin5s(:,k:k+1);
    TempMat=sortrows(TempMat,-1);
    for n=1:size(TempMat,2)
        h=imagesc(TempMat(:,n));set(gcf,'Visible', 'off');  
        colorbar
        caxis([min(TempMat(:)) max(TempMat(:))]);
        pngFileName = sprintf('Session_%d_%dbin_5s', (k+1)/2, n);
        saveas(gcf,fullfile('FiguresSpikes',pngFileName),format);
    end
end

%% Plotting PCo for the min by min analysis
%for 25ms
colormap jet
imagesc(CorrMat25);set(gcf,'Visible', 'off');
colorbar
saveas(gcf,fullfile('FiguresSpikes','PCo25ms'),'png');

%for 40ms
colormap jet
imagesc(CorrMat40);set(gcf,'Visible', 'off');
colorbar
saveas(gcf,fullfile('FiguresSpikes','PCo40ms'),'png');

%for 125ms
colormap jet
imagesc(CorrMat125);set(gcf,'Visible', 'off');
colorbar
saveas(gcf,fullfile('FiguresSpikes','PCo125ms'),'png');

%for 250ms
colormap jet
imagesc(CorrMat250);set(gcf,'Visible', 'off');
colorbar
saveas(gcf,fullfile('FiguresSpikes','PCo250ms'),'png');

%for 1s
colormap jet
imagesc(CorrMat1s);set(gcf,'Visible', 'off');
colorbar
saveas(gcf,fullfile('FiguresSpikes','PCo1s'),'png');

%for 5s
colormap jet
imagesc(CorrMat5s);set(gcf,'Visible', 'off');
colorbar
saveas(gcf,fullfile('FiguresSpikes','PCo5s'),'png');
toc;










