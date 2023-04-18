%% PAPER BDNF - PREDICTION ERROR (MMN) - MEG - BONETTI ET AL., HELIYON 2023

%%

%% PREPROCESSING - MEG DATA

%% MAXFILTER

%%% OBS!! REMEMBER TO CLOSE MATLAB.. THEN OPEN TERMINAL.. WRITE "use anaconda" AND THEN IN THE SAME TERMINAL OPEN MATLAB AGAIN AND SEND THE JOBS TO THE CLUSTER

maxfilter_path = '/neuro/bin/util/maxfilter';
project = 'MINDLAB2017_MEG-LearningBach';
maxDir = '/scratch7/MINDLAB2017_MEG-LearningBach/Portis/DIsp'; %output path

%actual construction of paths
path1 = '/raw/sorted/MINDLAB2015_MEG-TunteetMM';
path1d = dir([path1 '/0*']);
cnt = 2; %counter of mumufe people
for ii = 1:length(path1d) %over subject
    pd = dir([path1 '/' path1d(ii).name '/2*']);
    for dd = 1:length(pd)
        pd2 = dir([pd(dd).folder '/' pd(dd).name '/M*']);
        if strcmp(pd2.name,'MEG')
            pd3 = dir([pd2.folder '/' pd2.name '/0*']);
            for jj = 1:length(pd3)
                if strcmp(pd3(jj).name,'001.muMUFE_ns')
                    cnt = cnt + 1;
                    pd4 = dir([pd3(jj).folder '/' pd3(jj).name '/files/*fif']); %input fif file
                    rawName = [pd4.folder '/' pd4.name];
                    maxfName = [pd4.name(10:18) 'mumufe'];
                    %maxfilter to server.. (no movement compensation)
%                     cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,num2str(project), ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsss.fif'] ' -st 4 -corr 0.98 -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsss.log"']];
                    cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,num2str(project), ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsss.fif'] ' -st 4 -corr 0.98 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsss.log"']];
                    system(cmd);
                end
            end
        end
    end
end

%% addpath to OSL

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/ohba-external/fmt')
osl_startup

%% settings for cluster (parallel computing)

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
clusterconfig('slot', 6); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
% clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue

%% convertinig fif files from maxfilters into spm objects

fif_list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/DIsp/*fif');    % we made the list of sbj (fif_list) with the function dir

for ii = 1:length(fif_list)   %the loop will work for all ours sbj, so we put : from 1 to the length of fif_list
    S = [];                    % it's the structure
    S.dataset = [fif_list(ii).folder '/' fif_list(ii).name];  %we specified the name of the file that the loop has to take and we put (ii) to add one every time
%     jobid = job2cluster(@cluster_spmobject,S);
    D = spm_eeg_convert(S);     %it's our function for make the conversion. the structure S is he imput and D is the output 
end

%% REMOVING BAD TIME PERIODS USING OSLVIEW
%
% Now load oslview. This data has some exceptionally bad artefacts in. Mark the epochs at
% around 325s, 380s and 600s as bad, as well as everything from 650 seconds
% to the end. Marking is done by right-clicking in the proximity of the
% event and click on 'Mark Event'. A first click (green dashed label) marks
% the begin of a bad period, another second click indicates the end (in
% red). This will mean that we are not using about half of the data. But
% with such bad artefacts this is the best we can do. We can still obtain
% good results with what remains. NB: Push the disk button to save to disk
% (no prefix will be added, same name is kept).

%OBS!!! REMEMBER TO CHECK BAD SEGMENTS OF THE SIGNAL BOTH AT MEGPLANAR AND MEGMAG CHANNELS, YOU CAN CHANGE THE SELECTED CHANNELS THROUGH THE OSLVIEW INTERFACE

%OBS!!! REMEMBER TO MARK THE TRIAL WITHIN THE BAD SEGMENTS (BADTRIALS) AS
%BADTRIALS AND USE THE LABEL FOR REMOVING THEM FROM THE AVERAGING AFTER THE
%EPOCH!!

spm_list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/*mat');

for ii = 1:length(spm_list) % iterates over experimental blocks %OBS!!!!!
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    D = oslview(D);
    D.save(); %it is needed to save the marked bad events (and/or channels) in oslview
    disp(ii)
end

%% AFRICA denoising (ICA to remove eyeblink and heartbeat artefacts)

%changing channel labels..
%%% OBS!! to be ran only the first time %%%
for jj = 1:length(spm_list)
    D = spm_eeg_load([spm_list(jj).folder '/' spm_list(jj).name]);
    for ii = 1:length(D.chanlabels)
        idx = find(isspace(D.chanlabels{ii})==1); %removing space within data labels
        D = D.chanlabels(ii,[D.chanlabels{ii}(1:idx-1) D.chanlabels{ii}(idx+1:end)]); %making data labels equal to the ones of topoplot function for later ICA usage..
        disp([num2str(ii) '/' num2str(length(D.chanlabels)) ' - ' num2str(jj) '/' num2str(length(spm_list))])
    end
    D.save();
end
%%% until here %%%

%actual ICA calculation
% for ii = 77:length(spm_list) %:length(spm_files) %OBS!!!!! 
% v = [1,11,12,19,32,38,39,43];
% v = [38 39 43 50 70];
v = [11 12 19 32 39 43 50 70];
for ii = 1:length(v) 
    S = [];
    D = spm_eeg_load([spm_list(v(ii)).folder '/' spm_list(v(ii)).name]);
    S.D = D;
    jobid = job2cluster(@cluster_africa,S);
%     D = osl_africa(D,'do_ica',true,'do_ident',false,'do_remove',false,'used_maxfilter',true); 
%     D.save();
end

% v = [11 12 19 32];
v = [30];
%visual inspection and removal of artifacted components
for ii = 1:length(v)%:2:length(spm_list) %:length(spm_files) %OBS!!!!!
    D = spm_eeg_load([spm_list(v(ii)).folder '/' spm_list(v(ii)).name]);
    D = osl_africa(D,'do_ident',true,'do_remove',true);
    D.save();
    disp(v(ii))
end

%hacking the function to manage to get around the OUT OF MEMORY problem..
% jobid = job2cluster(@cluster_rembadcomp,S);

for ii = 31:length(spm_list)
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    D.montage('switch',1);
    disp(ii)
end

%D.montage('switch',1) %for seeing data before and after the AFRICA denoising

%% epoching (baseline = (-)100ms)

spm_list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/spm*');
prefix_tobeadded = 'e400'; %prefix to be added

%iterates for the recognition files only
for ii = 2:2:length(spm_list) %indexing only the files wanted for this paper
    %load the spm object
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    events = D.events;
    S = [];
    S.D = D;
    S.timewin = [-100 400];
    %event definitions
    S.trialdef(1).conditionlabel = 'Standard';
    S.trialdef(1).eventtype = 'STI 014_up';
    S.trialdef(1).eventvalue = 1:19;
    S.trialdef(2).conditionlabel = 'Pitch';
    S.trialdef(2).eventtype = 'STI 014_up';
    S.trialdef(2).eventvalue = 23:35;
    S.trialdef(3).conditionlabel = 'Timbre';
    S.trialdef(3).eventtype = 'STI 014_up';
    S.trialdef(3).eventvalue = 42:54;
    S.trialdef(4).conditionlabel = 'Localization';
    S.trialdef(4).eventtype = 'STI 014_up';
    S.trialdef(4).eventvalue = 61:73;
    S.trialdef(5).conditionlabel = 'Intensity';
    S.trialdef(5).eventtype = 'STI 014_up';
    S.trialdef(5).eventvalue = 80:92;
    S.trialdef(6).conditionlabel = 'Slide';
    S.trialdef(6).eventtype = 'STI 014_up';
    S.trialdef(6).eventvalue = 99:111;
    S.trialdef(7).conditionlabel = 'Rhythm';
    S.trialdef(7).eventtype = 'STI 014_up';
    S.trialdef(7).eventvalue = 118:130;
    S.prefix = prefix_tobeadded;
    De = spm_eeg_load([spm_list(ii).folder '/e400' spm_list(ii).name]);
    De.epochinfo = S;
    De.save();
    jobid = job2cluster(@cluster_epoch,S);
    %creates the epochinfo structure that is required for the source reconstruction later
%     epochinfo.trl = trl_sam;
%     epochinfo.time_continuous = D.time;
%     %switch the montage to 0 because for some reason OSL people prefer to do the epoching with the not denoised data
%     D = D.montage('switch',0);
end

%% assigning bad trials, if any..

spm_list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/spm*');
spm_list2 = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/e*');
for ii = 2%:2:length(spm_list) %indexing only the files wanted for this paper
    D2 = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]); %loading continuous data with proper D2.event field
    events = D2.events;
    %recreating proper event list..
    clear trl_sec
    ck = 0;
    for kkk = 1:length(events) %over events
        if strcmp(events(kkk).type,'STI 014_up') %with only STI 014_up
            ck = ck + 1;
            trl_sec(ck,1) = events(kkk).time;
            trl_sec(ck,2) = events(kkk).time + 0.3996;
        end
    end
    D = spm_eeg_load([spm_list2(ii).folder '/' spm_list2(ii).name]); %loading epoched data
    D = D.montage('switch',1);
    %take bad segments registered in oslview and check if they overlap with the trials. if so, it gives the number of overlapped trials that will be removed later   
    count = 0;
    Bad_trials = zeros(size(trl_sec,1),1);
    for kkk = 1:length(events) %over events
        if strcmp(events(kkk).type,'artefact_OSL')
            for k = 1:length(trl_sec) %over trials
                if events(kkk).time - trl_sec(k,2) < 0 %if end of trial is > than beginning of artifact
                    if trl_sec(k,1) < (events(kkk).time + events(kkk).duration) %if beginning of trial is < than end of artifact
                        Bad_trials(k,1) = 1; %it is a bad trial (stored here)
                        count = count + 1;
                    end
                end                  
            end
        end
    end
    %if bad trials were detected, their indices are stored within D.badtrials field
    if count == 0
        disp('there are no bad trials marked in oslview for');
        disp(ii/2);
    else
        D = badtrials(D,find(Bad_trials),1); %get the indices of the badtrials marked as '1' (that means bad)
        D.save(); %saving on disk
        disp('bad trials are ')
        length(D.badtrials)
    end
    epochinfo = [];
    trl_sam = round(trl_sec .* 600);
    trl_sam(:,2) = trl_sam(:,1) + 0.5*600;
    trl_sam(:,3) = -60;
    epochinfo.trl = trl_sam;
    epochinfo.time_continuous = D2.time;
    epochinfo.conditionlabels = D.conditions;
    D.epochinfo = epochinfo;
    D.save();
    D2.epochinfo = epochinfo;
    D2.save();
end

%% averaging and combining planar gradiometers (parallel computing)

list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/e*');

%settings for cluster (parallel computing)
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
clusterconfig('slot', 1); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
% clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue

%averaging
list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/esp*mat');
% addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Francesco/Functions'); %path to LEiDA_MEG_leonardo functions
output_prefix_to_be_set = 'm'; 
for ii = 2:2:length(list) %over files
    %distribute 
    input = [];
    input.D = [list(ii).folder '/' list(ii).name];
    input.prefix = output_prefix_to_be_set;
    jobid = job2cluster(@sensor_average, input); % this is the command for send the job to the cluster, in the brackets you can find the name on the function to run (afeter the @) and the variable for the input (in this case input)
    % look the script for more details about the function work
end
%combining planar gradiometers
list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/esp*mat');
for ii = 2:2:length(list) %over files
    input = [];
    input.D = [list(ii).folder '/' list(ii).name];
    D = spm_eeg_load(input.D);
    D = D.montage('switch',1);
    D.save();
    jobid = job2cluster(@combining_planar_cluster, input); % this is the command for send the job to the cluster, in the brackets you can find the name on the function to run (afeter the @) and the variable for the input (in this case input)
end

%% LBPD_startup_D

%starting up some of the functions that I wrote for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
addpath(pathl);
LBPD_startup_D(pathl);

%% extracting MEG sensor data

list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/Pmfh30*');

load_data = 1; %set 1 if you want to load the data instead of extracting it from SPM objects
datadir = '/scratch7/MINDLAB2017_MEG-LearningBach/Portis'; %path to data
outdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Portis/MEGsensor'; %path to write output in
% outdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper1/MCS_Sensor/lowpass40Hz/WM_RTs_Complexity'; %path to write output in
S = [];
%computing data
S.outdir = outdir;
S.data = [];
if load_data == 1 %if you already computed and saved on disk the t-tests you can load them here
    load([outdir '/MMNall_standardsubtracted.mat']);
    S.xlim = [0.22 0.24]; %time topolot
    S.topocondsing = [6];
    S.zlimmag = []; %magnetometers amplitude topoplot limits
    S.zlimgrad = [0 3]; %gradiometers amplitude topoplot limits
    S.data = data_mat;
    S.chanlabels = chanlabels;
    S.time_real = time_sel;
else %otherwise you can extract the data from SPM MEEG objects (one for each subject)
    S.spm_list = cell(1,length(list)/2);
    for ii = 2:2:length(list)
        S.spm_list(ii/2) = {[datadir '/' list(ii).name]};
    end
end
% S.conditions = {'Standard','Pitch','Slide','Intensity','Localization','Timbre','Rhythm'};
S.conditions = {'Pitch','Slide','Intensity','Localization','Timbre','Rhythm'};
% S.timeextract = [286:751]; %time-points to be extracted
S.timeextract = []; %time-points to be extracted
S.centerdata0 = 0; %1 to make data starting at 0
S.save_data = 1; %only meaningfull if you read data from SPM objects saved on disk
S.save_name_data = 'MMNall_standardsubtracted';
% S.save_name_data = 'sensor_data_lowpass3Hz_shortbaseline';
%individual waveform plotting
S.waveform_singlechannels_label = 0; %1 to plot single channel waveforms
%to plot waveform for deviants independently from one another
% S.condnumb = 6;
% S.condwaveform = {S.conditions{S.condnumb}};
S.color_line = [0.543 0 0];
S.wave_plot_conditions_together = 0; %1 for plotting the average of all conditions
S.mag_lab = 2; %1 for magnetometers; 2 for gradiometers
S.x_lim_temp_wave = []; %limits for time (in secs) (E.g. [-0.1 3.4])
S.y_lim_ampl_wave = [-1 4]; %limit for amplitude (E.g. [0 120] magnetometes, [0 6] gradiometers)
%averaged waveform plotting
S.waveform_average_label = 0;
% S.left_mag = [2:2:204];
S.legc = 1; %set 1 for legend
S.left_mag = 98;
S.signtp = {[]};
S.avewave_contrast = 0; %1 to plot the contrast between conditions (averaged waveform)
S.save_label_waveaverage = 0;
S.label_plot = 'MMN';
%t-tests
S.t_test_for_permutations = 0;
S.cond_ttests_tobeplotted_topoplot = [1 2]; %this is for both topoplot and t-tests!! (here [1 2] means cond1 vs cond2!!!!!!!)
%topoplotting
S.topoplot_label = 1;
S.fieldtrip_mask = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External';
S.topocontr = 0;
S.colormap_spec = 0;
% x = []; x.bottom = [0 0 1]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 1 0.5]; x.top = [1 0.95 0]; %yellow - blue
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
S.colormap_spec_x = x;
S.topoplot_save_label = 0;

[out] = MEG_sensors_plotting_ttest_LBPD_D(S);

%% subtracting standard from MMN (deviants) and saving it..

JK = zeros(size(data_mat,1),size(data_mat,2),size(data_mat,3),size(data_mat,4)-1);
cnt = 0;
for ii = 2:7
    cnt = cnt + 1;
    JK(:,:,:,cnt) = data_mat(:,:,:,ii)-data_mat(:,:,:,1);
end
data_mat = JK;
save MMNall_standardsubtracted_f30.mat data_mat S time_sel chanlabels

%%

%% (BETTER) ANALYSES PERFORMED FOR THE REVISION IN HELYION

%extracting peaks of MMN for maximum amplitude MEG channels, then computing
%2-way ANOVAs (musicianship and BDNF as factors) and MMN amplitude as dependent variable 

%alway run this section first
%then run the section below that you need

load('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/Criscuolo2019_mus_vs_nonm.mat') %loading groups
load('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/MEGsensor/MMNall_standardsubtracted.mat'); %loading data
load('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/groups.mat')
%peak channels (left): 9 11 15 109 (211 221 241 1511); (right): 93 95 97 197 (1311 1321 1331 2611)
%peaks: -pit: 170-190 (163-175); sli: 150-170 (151-163); int: 130-150 (139-151); loc: 110-130 (127-139); tim: 125-145 (136-148); rhy: 260-280 (217-229)
% S.conditions = {'Pitch','Slide','Intensity','Localization','Timbre','Rhythm'};
data_mat2 = abs(data_mat(:,:,A,:)); %getting only participants with genetic data (abs because we have double polarity in magnetometers)
%new subject IDs divided into groups
clear NM1 NM0 M1 M0
cnm1 = 0; cnm0 = 0; cm1 = 0; cm0 = 0;
for ii = 1:length(char_group_idx_fin)
    if strcmp(char_group_idx_fin{ii},'nonm1')
        cnm1 = cnm1 + 1;
        NM1(cnm1) = ii;
    elseif strcmp(char_group_idx_fin{ii},'nonm0')
        cnm0 = cnm0 + 1;
        NM0(cnm0) = ii;
    elseif strcmp(char_group_idx_fin{ii},'mus1')
        cm1 = cm1 + 1;
        M1(cm1) = ii;
    else
        cm0 = cm0 + 1;
        M0(cm0) = ii;
    end
end
data_left = zeros(length(char_group_idx_fin),6);
data_right = zeros(length(char_group_idx_fin),6);
data_tot = zeros(length(char_group_idx_fin),6);
timew{1} = 163:175; timew{2} = 151:163; timew{3} = 139:151; timew{4} = 127:139; timew{5} = 136:148; timew{6} = 217:229;
for ii = 1:6 %over deviants
    data_left(:,ii) = squeeze(mean(mean(data_mat2([9 11 15 109],timew{ii},:,ii),2),1)); %data averaged over peaks and peak-channels for every participant (74) and deviant (6) (left hemisphere
    data_right(:,ii) = squeeze(mean(mean(data_mat2([93 95 97 197],timew{ii},:,ii),2),1)); %same but for right hemisphere
    data_tot(:,ii) = squeeze(mean(mean(data_mat2([9 11 15 109 93 95 97 197],timew{ii},:,ii),2),1)); %magnetometers same but for right hemisphere
end

%% two-way ANOVA (anovan Matlab function)

dev = 6; % {'Pitch','Slide','Intensity','Localization','Timbre','Rhythm'};
leftl = 2; %1 = left; 0 = right; 2 = left and right averaged together

if leftl == 1
    dataano = data_left(:,dev);
elseif leftl == 0
    dataano = data_right(:,dev);
elseif leftl == 2
    dataano = data_tot(:,dev);
end

[p,t,stats] = anovan(dataano,{g1,g2},'model','interaction','varnames',{'g1','g2'}); %'off' for not showing the plot
% [c,m,h,nms] = multcompare(stats,'ctype','tukey-kramer'); %perform multiple comparison test based on anova1 output = post-hoc analysis (c needs to be saved for every time point and every channel)
clc
if leftl == 1
    disp('M1')
    mean(data_left(M1,dev),1)
    disp('M0')
    mean(data_left(M0,dev),1)
    disp('NM1')
    mean(data_left(NM1,dev),1)
    disp('NM0')
    mean(data_left(NM0,dev),1)
elseif leftl == 0
    disp('M1')
    mean(data_right(M1,dev),1)
    disp('M0')
    mean(data_right(M0,dev),1)
    disp('NM1')
    mean(data_right(NM1,dev),1)
    disp('NM0')
    mean(data_right(NM0,dev),1)
elseif leftl == 2
    disp('M1')
    mean(data_tot(M1,dev),1)
    std(data_tot(M1,dev),[],1)
    disp('M0')
    mean(data_tot(M0,dev),1)
    std(data_tot(M0,dev),[],1)
    disp('NM1')
    mean(data_tot(NM1,dev),1)
    std(data_tot(NM1,dev),[],1)
    disp('NM0')
    mean(data_tot(NM0,dev),1)
    std(data_tot(NM0,dev),[],1)
end

%plotting reaction times
%some default color specifications for later plotting..
% cb = [0.5 0.8 0.9; 1 1 0.7; 0.7 0.8 0.9; 0.8 0.5 0.4; 0.5 0.7 0.8; 1 0.8 0.5; 0.7 1 0.4; 1 0.7 1; 0.6 0.6 0.6; 0.7 0.5 0.7; 0.8 0.9 0.8; 1 1 0.4];
cl = [0.7 0 0; 0.7 0 0; 0 0 0.8; 0 0 0.8; 0 0 0.8;];
datadum = cell(1,4);
datadum{1} = data_tot(M1,dev);
datadum{2} = data_tot(M0,dev);
datadum{3} = data_tot(NM1,dev);
datadum{4} = data_tot(NM0,dev);
figure
rm_raincloud2(datadum',cl)
grid minor
set(gcf,'color','w')
set(gcf,'Position',[200,200,400,550])
xlabel('Amplitude (fT)'); %set(gca,'YDir','normal');
% xlim([1500 3500])

%% Correlations (not reported in the paper but asked during the revision)

dev = 1; % {'Pitch','Slide','Intensity','Localization','Timbre','Rhythm'};
groupl = 4; %1 = M1; 2 = M0; 3 = NM1; 4 = NM0 

clc
load('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/MEGsensor/Stats/musicaltraining.mat');
leftl = 2; %1 = left; 0 = right; 2 = left and right averaged together
if leftl == 1
    dataano = data_left(:,dev);
elseif leftl == 0
    dataano = data_right(:,dev);
elseif leftl == 2
    dataano = data_tot(:,dev);
end
if groupl == 1
    ghj = M1;
elseif groupl == 2
    ghj = M0;
elseif groupl == 3
    ghj = NM1;
elseif groupl == 4
    ghj = NM0;
end

[rho,p] = corr(dataano(ghj),azz(ghj)) % correlation function 
% [rho,p] = corr(dataano,azz) % correlation function

%% testing significance of difference in musical training in the four groups

load('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/groups.mat')
load('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/MEGsensor/Stats/musicaltraining.mat');
[p,t,stats] = anovan(azz,{g1,g2},'model','interaction','varnames',{'g1','g2'}); %'off' for not showing the plot

%% THE ADDITIONAL/NOVEL ANALYSES REQUESTED BY HELIYON"S REVIEWERS WHICH INTEGRATED/REPLACED THE ORIGINAL ONES ENDED HERE

%%

%% CONTRASTS IN SOURCES (JUST A REMINDER.. NOT REALLY USEFUL HERE..)

% oat.source_recon.conditions        = {'Standard','Pitch','Timbre','Localization','Intensity','Slide','Rhythm'};
% oat.first_level.contrast{1} = [1 0 0 0 0 0 0]'; %Standard
% oat.first_level.contrast{2} = [0 1 0 0 0 0 0]'; %Pitch
% oat.first_level.contrast{3} = [0 0 1 0 0 0 0]'; %Timbre
% oat.first_level.contrast{4} = [0 0 0 1 0 0 0]'; %Localization
% oat.first_level.contrast{5} = [0 0 0 0 1 0 0]'; %Intensity
% oat.first_level.contrast{6} = [0 0 0 0 0 1 0]'; %Slide
% oat.first_level.contrast{7} = [0 0 0 0 0 0 1]'; %Rhythm
% oat.first_level.contrast{8} = [-1 1 0 0 0 0 0]'; %Pitch - Standard  peak 173
% oat.first_level.contrast{9} = [-1 0 1 0 0 0 0]'; %Timbre - Standard peak 182
% oat.first_level.contrast{10} = [-1 0 0 1 0 0 0]'; %Localization - Standard peak 128
% oat.first_level.contrast{11} = [-1 0 0 0 1 0 0]'; %Intensity - Standard peak 152
% oat.first_level.contrast{12} = [-1 0 0 0 0 1 0]'; %Slide - Standard peak 158
% oat.first_level.contrast{13} = [-1 0 0 0 0 0 1]'; %Rhythm - Standard peak 158

%% GETTING ORIGINAL SUBJECT ID (useful for checking and calculating demographic data, etc.)

%loading information for getting proper ID of musicinas and BDNF
load('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/Criscuolo2019_mus_vs_nonm.mat') %load groups
load('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/S_for_Structure.mat'); %load previously saved settings
group = {mus0 nonm0 mus1 nonm1};
%getting proper IDs and having them as double
MUS0 = S.spm_list(group{1});
MUS1 = S.spm_list(group{3});
nMUS0 = S.spm_list(group{2});
nMUS1 = S.spm_list(group{4});
mus0b = zeros(length(MUS0),1);
for ii = 1:length(MUS0)
   mus0b(ii,1) = str2double(MUS0{ii}(64:67));   
end
mus1b = zeros(length(MUS1),1);
for ii = 1:length(MUS1)
   mus1b(ii,1) = str2double(MUS1{ii}(64:67));   
end
nmus0b = zeros(length(nMUS0),1);
for ii = 1:length(nMUS0)
   nmus0b(ii,1) = str2double(nMUS0{ii}(64:67));   
end
nmus1b = zeros(length(nMUS1),1);
for ii = 1:length(nMUS1)
   nmus1b(ii,1) = str2double(nMUS1{ii}(64:67));   
end

%% getting actual demographica data (TO BE RUN AFTER RUNNING THE SECTION ABOVE)

groupn = 2; %number of the group (1 = M1; 2 = M0; 3 = NM1; 4 = NM0; 5 = all subjects)
playn_train = 4; % 4 = formal musical training; 5 = musical playing

clc
group{1} = mus1b; group{2} = mus0b; group{3} = nmus1b; group{4} = nmus0b; group{5} = cat(1,mus1b,mus0b,nmus1b,nmus0b);
[~,~,raw] = xlsread('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/BDNF_MMN/DemographicalData.xlsx');
clear ASD
cnt = 0;
for ii = 1:140 %over subjects
    if find(str2double(raw{ii+3,1}(4:end)) == group{groupn}) %if the subject matches the group
        cnt = cnt + 1;
        ASD(cnt,1:5) = raw(ii+3,3:7); %storing subjects of the chosen group
    end
end
%mean age
ave = mean(cell2mat(ASD(:,1)));
std2 = std(cell2mat(ASD(:,1)));
disp(['total number of subjects is : ' num2str(length(ASD))])
disp(['mean age is : ' num2str(ave) ' / std is : ' num2str(std2)])
%mean years of musical training/playing (5 = formal musical training; 6 = musical playing)
ave = mean(cell2mat(ASD(:,playn_train)));
std2 = std(cell2mat(ASD(:,playn_train)));
disp(['mean age of musical training is : ' num2str(ave) ' / std is : ' num2str(std2)])
%handedness and sex
cL = 0;
cR = 0;
cA = 0;
cM = 0;
cF = 0;
for ii = 1:length(ASD) %over subjects of chosen group
    %handedness
    if strcmp(ASD(ii,2),'L')
        cL = cL + 1;
    elseif strcmp(ASD(ii,2),'R')
        cR = cR + 1;
    elseif strcmp(ASD(ii,2),'A')
        cA = cA + 1;
    end
    %sex
    if strcmp(ASD(ii,3),'M')
        cM = cM + 1;
    elseif strcmp(ASD(ii,3),'F')
        cF = cF + 1;
    end
end
disp(['males are : ' num2str(cM) ' / females are : ' num2str(cF)])
disp(['left-handed are : ' num2str(cL) ' / right-handed are : ' num2str(cR) ' / ambidextriuos are : ' num2str(cA)])

%%

%% SOURCE RECONSTRUCTION

%% DICOM to NIFTI conversion 

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/dicm2nii'); %adds path to the dcm2nii folder in osl
MRIoutput = '/scratch7/MINDLAB2017_MEG-LearningBach/Portis/MRI';
% MRIpath1 = '/raw/sorted/MINDLAB2015_MEG-TunteetMM';
MRIsubj = dir('/raw/sorted/MINDLAB2015_MEG-TunteetMM/0*');
cnt = 1;

for ii = 2:length(MRIsubj)
    %building path to dicom files
    MRIMEGdate = dir([MRIsubj(ii).folder '/' MRIsubj(ii).name '/20*']);
    for jj = 1:length(MRIMEGdate)
        MEGMRI = dir([MRIMEGdate(jj).folder '/' MRIMEGdate(jj).name '/M*']);
        if strcmp(MEGMRI.name(1:2),'MR')
            MEGMRI2 = dir([MEGMRI.folder '/' MEGMRI.name '/*ant']);
            if isempty(MEGMRI2)
                warning(['subj ' num2str(ii) ' is not ending with MPRAGE-variant..'])
            elseif length(MEGMRI2) > 1
                warning(['subj ' num2str(ii) ' has more than one MR folder..'])
            end
            dcmSource = [MEGMRI2(1).folder '/' MEGMRI2(1).name '/files/'];
            niiFolder = [MRIoutput '/' MRIsubj(ii).name];
            dicm2nii(dcmSource, niiFolder, '.nii'); %actual conversion from dicom to niftii
        end
    end
    cnt = cnt + 1;
    disp(ii)
end

%% settings for cluster (parallel computing)

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
clusterconfig('slot', 3); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
% clusterconfig('scheduler', 'none'); % set automatically the long run queue
clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue

%% RHINO COREGISTRATION (parallel computing)

a = dir(['/scratch7/MINDLAB2017_MEG-LearningBach/Portis/MRI/0*']);
spm_list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/spm*');

for ii = 2:2:length(spm_list)
    S = [];
    S.ii = ii;
    S.D1 = [spm_list(ii).folder '/' spm_list(ii).name];
    S.D = [spm_list(ii).folder '/e400' spm_list(ii).name]; %check this if you want 'epoched' or 'continuous'
    cnt2 = 1;
    while ~strcmp(spm_list(ii).name(12:15),a(cnt2).name) && cnt2 < length(a)
        cnt2 = cnt2 + 1;
    end
    if strcmp(spm_list(ii).name(12:15),a(cnt2).name) || cnt2 ~= 119
        caz = dir([a(cnt2).folder '/' a(cnt2).name '/M*']);
        S.mri = [caz(1).folder '/' caz(1).name];
    else
        S.mri = ['/scratch7/MINDLAB2017_MEG-LearningBach/Portis/' spm_list(ii).name(12:15) 'MRI_MNI152/MNI152_T1_2mm.nii'];
    end
    S.useheadshape = 1;
    S.use_rhino = 1; %set 1 for having rhino, 0 for not having it
    S.forward_meg = 'Single Shell';
    S.fid.label.nasion = 'Nasion';
    S.fid.label.lpa = 'LPA';
    S.fid.label.rpa = 'RPA';

    jobid = job2cluster(@coregfuncp,S); %running with parallel computing
    
    disp(num2str(ii))
end

%% checking results of rhino coregistration

D = spm_eeg_load(S.D);
D1 = spm_eeg_load(S.D1);
D = D.montage('switch',1); %switch the montage to 1 in order to be safer
rhino_display(D);
D1 = D1.montage('switch',1); %switch the montage to 1 in order to be safer
rhino_display(D1);


%% back to before AFRICA (for MEG channels label issues..)
     
D2 = spm_eeg_load('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/DIsp/spmeeg_SUBJ0002_mumufe_tsss.mat');

spm_list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/es*');
v = [64 70 136 142 146 154 230];
for ii = 1:length(v)%60%4:2:length(spm_list)
    D = spm_eeg_load([spm_list(v(ii)).folder '/' spm_list(v(ii)).name]);
    S2 = [];
    S2.D = D; %file to be updated progressively
    S2.D2 = D2; %original labels
    jobid = job2cluster(@cluster_Dlabel,S2); %running with parallel
end

%% settings for cluster (parallel computing)

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
clusterconfig('slot', 6); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
% clusterconfig('scheduler', 'none'); % set automatically the long run queue
clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue
    
%% ACTUAL CODES FOR SOURCE RECONSTRUCTION - BEAMFORMING (OAT CODES)

spm_list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/es*');

v = [64 70 136 142 146 154 230];
for ii = 1:length(v) %6:2:length(spm_list)
    processed_file = [spm_list(v(ii)).folder '/' spm_list(v(ii)).name];
    D_epoched = spm_eeg_load(processed_file);
    D_epoched = D_epoched.montage('switch',1); %switch the montage to 1 in order to be safer (so we have the AFRICA denoised data)
    D_epoched.save();
    % Beamform
    oat = [];
    oat.source_recon.D_epoched(str2double(spm_list(v(ii)).name(13:16)))         = {processed_file};
%     oat.source_recon.D_epoched(str2double(spm_list(v(ii)).name(14:17)))         = {processed_file};
    pca_dim_1 = 50;
    oat.source_recon.pca_dim = pca_dim_1(ones(length(oat.source_recon.D_epoched),1)); %this seems to be necessary. It is for having the same number of values than the number of input D objects
    oat.source_recon.modalities = {'MEGMAG'; 'MEGPLANAR'};
    oat.source_recon.conditions        = {'Standard','Pitch','Timbre','Localization','Intensity','Slide','Rhythm'};
    oat.source_recon.gridstep          = 8; % in mm
    oat.source_recon.time_range        = [-0.1 0.39]; % time range in secs
    oat.source_recon.freq_range        = [0.1 40]; % frequency range in Hz
%     oat.source_recon.freq_range        = [0.1 10]; % frequency range in Hz
    %S.source_recon.pca_order         = 250;
    oat.source_recon.type              = 'Scalar';
    oat.source_recon.method            = 'beamform';
    oat.source_recon.normalise_method  = 'mean_eig';
    oat.source_recon.forward_meg       = 'Single Shell';
    %S.source_recon.prefix            = '';
    oat.source_recon.report.do_source_variance_maps = 1;
    oat.source_recon.sessions_to_do = [];
    oat.source_recon.sessions_to_do = [str2double(spm_list(v(ii)).name(13:16))]; %sessions to do among the file_list (subject 39 excluded because of problems during data collection..)
%     oat.source_recon.sessions_to_do = [str2double(spm_list(v(ii)).name(14:17))]; %sessions to do among the file_list (subject 39 excluded because of problems during data collection..)
    oat.source_recon.dirname           = [spm_list(1).folder '/sourcetryindsubj600Hz/' spm_list(v(ii)).name(13:16)];
%     oat.source_recon.dirname           = [spm_list(1).folder '/sourcetryindsubj600Hz/' spm_list(v(ii)).name(14:18)];   
    jobid = job2cluster(@cluster_beamforming,oat); %running with parallel
end

%% PROCEEDING WITH FIRST LEVEL ANALYSIS, ETC. OSL - OAT

%% codes to get a "normal" oat for source reconstruction (so all subjects in the same folder)..
%this is a trick to be able to use the parallel computing of our cluster

spm_list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/es*');

% v = [64 70 136 142 146 154 230];
oat = [];
for ii = 2:2:length(spm_list)
    processed_file = [spm_list(ii).folder '/' spm_list(ii).name];
    oat.source_recon.D_epoched(ii/2)         = {processed_file};
    oat.source_recon.sessions_to_do(ii/2) = {[str2double(spm_list(ii).name(13:16))]}; %sessions to do among the file_list (subject 39 excluded because of problems during data collection..)
    oat.source_recon.results_fnames(ii/2) = {['session' num2str(str2double(spm_list(ii).name(13:16))) '_recon']};
end
% Beamform
% oat.source_recon.D_epoched(str2double(spm_list(ii).name(13:16)))         = {processed_file};
pca_dim_1 = 50;
oat.source_recon.pca_dim = pca_dim_1(ones(length(oat.source_recon.D_epoched),1)); %this seems to be necessary. It is for having the same number of values than the number of input D objects
oat.source_recon.modalities = {'MEGMAG'; 'MEGPLANAR'};
oat.source_recon.conditions        = {'Standard','Pitch','Timbre','Localization','Intensity','Slide','Rhythm'};
oat.source_recon.gridstep          = 8; % in mm
oat.source_recon.time_range        = [-0.1 0.39]; % time range in secs
oat.source_recon.freq_range        = [0.1 40]; % frequency range in Hz
oat.source_recon.type              = 'Scalar';
oat.source_recon.method            = 'beamform';
oat.source_recon.normalise_method  = 'mean_eig';
oat.source_recon.forward_meg       = 'Single Shell';
oat.source_recon.report.do_source_variance_maps = 1;
oat.source_recon.dirname           = [spm_list(1).folder '/sourcetryindsubj600Hz/firstlevel'];

%% moving files computed independently for each subject into a common new folder (trying to go back to the "usual" way of running oat analysis pipeline in OSL)..

%%% TO BE RUN ONLY THE FIRST TIME %%%
list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/0*');
bs = '/scratch7/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/firstlevel.oat';
for ii = 1:length(list)
    %concat file .dat
    as = ['/scratch7/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/' list(ii).name '/concatMfsession' num2str(str2double(list(ii).name(1:4))) '_spm_meeg.dat'];
    status = movefile(as,bs)
    %concat file .mat
    as = ['/scratch7/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/' list(ii).name '/concatMfsession' num2str(str2double(list(ii).name(1:4))) '_spm_meeg.mat'];
    status = movefile(as,bs)
    %session file
    as = ['/scratch7/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/' list(ii).name '/session' num2str(str2double(list(ii).name(1:4))) '_recon.mat'];
    status = movefile(as,bs)
end

%% settings for cluster (parallel computing)

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
clusterconfig('slot', 1); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
% clusterconfig('scheduler', 'none'); % set automatically the long run queue
clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue

%% first level (each experimental block for each subject, independently)

%FIRST LEVEL
design_matrix_summary = {};
design_matrix_summary{1} = [1 0 0 0 0 0 0];design_matrix_summary{2} = [0 1 0 0 0 0 0]; design_matrix_summary{3}=[0 0 1 0 0 0 0];design_matrix_summary{4}=[0 0 0 1 0 0 0]; design_matrix_summary{5}=[0 0 0 0 1 0 0]; design_matrix_summary{6}=[0 0 0 0 0 1 0]; design_matrix_summary{7}=[0 0 0 0 0 0 1];
oat.first_level.design_matrix_summary = design_matrix_summary;
% contrasts to be calculated:
oat.first_level.contrast = {};
% contrast design matrix
oat.source_recon.conditions        = {'Standard','Pitch','Timbre','Localization','Intensity','Slide','Rhythm'};
oat.first_level.contrast{1} = [1 0 0 0 0 0 0]'; %Standard
oat.first_level.contrast{2} = [0 1 0 0 0 0 0]'; %Pitch
oat.first_level.contrast{3} = [0 0 1 0 0 0 0]'; %Timbre
oat.first_level.contrast{4} = [0 0 0 1 0 0 0]'; %Localization
oat.first_level.contrast{5} = [0 0 0 0 1 0 0]'; %Intensity
oat.first_level.contrast{6} = [0 0 0 0 0 1 0]'; %Slide
oat.first_level.contrast{7} = [0 0 0 0 0 0 1]'; %Rhythm
oat.first_level.contrast{8} = [-1 1 0 0 0 0 0]'; %Pitch - Standard
oat.first_level.contrast{9} = [-1 0 1 0 0 0 0]'; %Timbre - Standard
oat.first_level.contrast{10} = [-1 0 0 1 0 0 0]'; %Localization - Standard
oat.first_level.contrast{11} = [-1 0 0 0 1 0 0]'; %Intensity - Standard
oat.first_level.contrast{12} = [-1 0 0 0 0 1 0]'; %Slide - Standard
oat.first_level.contrast{13} = [-1 0 0 0 0 0 1]'; %Rhythm - Standard
%contrast names
oat.first_level.contrast_name = {};
oat.first_level.contrast_name{1} = 'Standard';
oat.first_level.contrast_name{2} = 'Pitch';
oat.first_level.contrast_name{3} = 'Timbre';
oat.first_level.contrast_name{4} = 'Localization';
oat.first_level.contrast_name{5} = 'Intensity';
oat.first_level.contrast_name{6} = 'Slide';
oat.first_level.contrast_name{7} = 'Rhythm';
oat.first_level.contrast_name{8} = 'Pitch - Standard';
oat.first_level.contrast_name{9} = 'Timbre - Standard';
oat.first_level.contrast_name{10} = '%Localization - Standard';
oat.first_level.contrast_name{11} = 'Intensity - Standard';
oat.first_level.contrast_name{12} = 'Slide - Standard';
oat.first_level.contrast_name{13} = 'Rhythm - Standard';
%original contrasts, try1 and try2 seem to give the exact same results
oat.first_level.report.first_level_cons_to_do = [1 2 3 4 5 6 7 8 9 10 11 12 13]; %[1 2 3]; %better to do 3 2 1 in order to get the information for the peak value of the contrast old-new
oat.first_level.time_range = [-0.1 0.39];
oat.first_level.post_tf_downsample_factor = 1;
%slow negativity
% oat.first_level.cope_type = 'coape'; %this should calculate absolute values of glm copes before comparing different conditions (otherwise, by default the absolute value is done on the copes after comparing different conditions); the absolute value here seems to be necessary because of the sign ambiguity issue introduced by source reconstruction
%N100
oat.first_level.cope_type = 'none'; %this should calculate absolute values of glm copes before comparing different conditions (otherwise, by default the absolute value is done on the copes after comparing different conditions); the absolute value here seems to be necessary because of the sign ambiguity issue introduced by source reconstruction
oat.first_level.name = ['wholebrain_first_level_BC_nocoape']; %REMEMBER TO CHECK THIS NAME!!
oat.first_level.bc = ones(1,13);
%to add if the oat has not been automatically saved
list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/0*');
for ii = 1:length(list)
%     oat.first_level.results_fnames(ii) = {['session' num2str(str2double(list(ii).name(1:4))) '_wholebrain_first_level_BC.mat']};
    oat.first_level.results_fnames(ii) = {['session' num2str(str2double(list(ii).name(1:4))) '_wholebrain_first_level_BC_nocoape.mat']};
end

%% running first level on parallel computing
for ii = 1:length(list)
    oat.first_level.sessions_to_do = [];
    oat.first_level.sessions_to_do = [ii]; %here it seems that the session indexes oat.source_recon.results_fnames{ii} is directly related to the sequential 
%     oat.first_level.sessions_to_do = [str2double(list(ii).name(1:4))];
    jobid = job2cluster(@cluster_beamfirstlevel,oat);
end

%% SUBJECT LEVEL
%in this case it does not do anything since I have only one experimental
%block for each participant (however, for computational reasons, I have to
%run it)
%this is needed to read the proper subjects..
for ii = 1:140
%     oat.first_level.results_fnames(ii) = {['session' num2str(ii) '_wholebrain_first_level_BC.mat']};
    oat.first_level.results_fnames(ii) = {['session' num2str(ii) '_wholebrain_first_level_BC_nocoape.mat']};
end
% this is if you have a perfect correspondance between sessions and subjects
oat.subject_level.session_index_list = cell(1,140); %sarebbe length(subjects)
oat.subject_level.name = 'MMN';
oat.subject_level.subjects_to_do = [];
% oat.subject_level.subjects_to_do = [5:115];
for ii = 1:length(list)
    oat.subject_level.subjects_to_do(ii) = str2double(list(ii).name(1:4));
end
%to update the names of the results files of subject level
for ii = 1:140
    oat.subject_level.session_index_list{ii} = ii;
    oat.subject_level.results_fnames(ii) = {['subject' num2str(ii) '_wholebrain_first_level_BC_nocoape_MMN.mat']};
end

oat = osl_check_oat(oat);
oat.to_do = [0 0 1 0];
oat = osl_run_oat(oat);

%% GROUP LEVEL - MAIN EFFECT OF MMN SOURCE RECONSTRUCTION

oat.group_level = [];
oat.group_level.name = 'group_level_everybody_BC'; %OBS!! REMEMBER TO UPDATE THE NAME!
oat.group_level.subjects_to_do = [];
% oat.group_level.subjects_to_do = 1:70; %this seems to be necessary since the function gets the indices for the group level from each progressive number of oat.subject_level.subjects_to_do (not very clever idea by the way..). Therefore not to mix up things, you should have here a vector from 1 to the total number of subjects and then specify the correct subject IDs that you want.. or at least, to me it looks like that..
oat.group_level.subjects_to_do = 1:115; %this seems to be necessary since the function gets the indices for the group level from each progressive number of oat.subject_level.subjects_to_do (not very clever idea by the way..). Therefore not to mix up things, you should have here a vector from 1 to the total number of subjects and then specify the correct subject IDs that you want.. or at least, to me it looks like that..
%results name
% oat.group_level.results_fnames = ['wholebrain_first_level_BC_MMN' '_' oat.group_level.name '.mat'];
oat.group_level.results_fnames = ['wholebrain_first_level_BC_nocoape_MMN' '_' oat.group_level.name '.mat'];
% Spatial and temporal averaging options
oat.group_level.time_range = [-0.1 0.39];
oat.group_level.space_average = 0;
oat.group_level.time_average = 0;
oat.group_level.time_smooth_std = 0; % secs
oat.group_level.use_tstat = 0;
%path to AAL template
% oat.group_level.mask_fname = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_8mm_try5.nii.gz';
% Spatial and temporal smoothing options
oat.group_level.spatial_smooth_fwhm = 0; % mm
oat.group_level.group_varcope_time_smooth_std = 100;
oat.group_level.group_varcope_spatial_smooth_fwhm = 100; % smooths the variance of the group copes. It is recommended to do this.
%store copes (useful for doing the permutation test later, otherwise it needs to compute again the group level analysis)
oat.group_level.store_lower_level_copes = 1;
% Set up design matrix and contrasts
%this is if you have only the general mean across the all participants
oat.group_level.group_design_matrix = ones(1,length(oat.group_level.subjects_to_do)); %if you want all of the participants
oat.group_level.group_contrast = [];
oat.group_level.group_contrast{1} = [1];
oat.group_level.group_contrast_name = {};
oat.group_level.group_contrast_name{1} = 'mean';
oat.group_level.glm_method='fixed_effects'; %ols or fixed-effects
% Define which contrasts to perform for the report
oat.group_level.first_level_contrasts_to_do = [1,2,3,4,5,6,7,8,9,10,11,12,13]; % list of first level contrasts to run the group analysis on
oat.group_level.report.first_level_cons_to_do = [1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.group_level_cons_to_do = [1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.show_lower_level_copes = 0;
oat.group_level.report.show_lower_level_cope_maps = 0;

%% creatig nifti images with statistics

load('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/firstlevel.oat/oat_wholebrain_first_level_BC_MMN_group_level_everybody_BC.mat');

S2 = [];
S2.oat = oat;
S2.stats_fname = oat.group_level.results_fnames;
S2.first_level_contrasts = 1:13; %remember that in this way you define the order of the output (this numbers refers to the order (numbers) defined in the contrasts; for example here tstat3 refers to the contrast old-new)
S2.group_level_contrasts = [1 2];
S2.resamp_gridstep = oat.source_recon.gridstep;

jobid = job2cluster(@cluster_oat_save_nii_stats,S2);

%% settings for cluster (parallel computing)

addpath('/projects/MINDLAB2020_MEG-AuditoryPatternRecognition/scripts/leonardo')
clusterconfig('slot', 2); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
% clusterconfig('scheduler', 'none');
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); % 0 = short; 1 = all; 2 = long


%% clustering in source space (peak MMN)

%PEAKS for the deviants (ms from onset)
%int = 152; loc = 128; pit = 173; rhy = 158; sli = 158; tim = 182

%preparing settings
%loading proper oat (mean over subjects)
load('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/firstlevel.oat/oat_wholebrain_first_level_BC_MMN_group_level_everybody_BC.mat');
outdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/BDNF_MMN';
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster') %add the path where is the function for submit the jobs to the server
% load('/scratch7/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/firstlevel.oat/oat_wholebrain_first_level_BC_MMN_group_level_mus1vsmus0_BC.mat');
dev = 8:13; %deviants
timep = cell(1,6);
timep{4} = [0.142 0.162]; %time-range for cluster calculation (in seconds)
timep{3} = [0.118 0.138]; %time-range for cluster calculation (in seconds)
timep{1} = [0.163 0.183]; %time-range for cluster calculation (in seconds)
timep{6} = [0.148 0.168]; %time-range for cluster calculation (in seconds)
timep{5} = [0.148 0.168]; %time-range for cluster calculation (in seconds)
timep{2} = [0.172 0.192]; %time-range for cluster calculation (in seconds)
t_val = 2; %t-value to binarize statistics (in the paper we indicated )

%actual computation
for ii = 2:6
    S = [];
    % oat.clustname = ['tone_' num2str(ii)];
    S.oat = oat; %oat computed in this section
    S.cluster_stats_thresh = t_val; %t-value to binarize statistics
    S.cluster_stats_nperms = 5000; %permutations
    S.first_level_copes_to_do = oat.group_level.first_level_contrasts_to_do(dev(ii));
    S.group_level_copes_to_do = 1; %MUS1 > MUS0
    S.group_varcope_spatial_smooth_fwhm = oat.group_level.group_varcope_spatial_smooth_fwhm;
    S.write_cluster_script = 0;
    S.time_range = timep{ii};
    S.time_average=1;
    % Run the permutations (on cluster.. parallel computing)
    jobid = job2cluster(@clusterbasedpermutation_osl,S);
end

%% GETTING INFORMATION FROM THE SIGNIFICANT SOURCE CLUSTERS (AND FOR THE GENERAL MMN SOURCES)

%creating a file with information on significant clusters/voxels outputted by the permutation test (THINK TO MAKE THIS A FUNCTION)
out = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/BDNF_MMN';

%actual name plus path
fg = dir([out '/*generalMMN_source.nii.gz']);
for mm = 1:length(fg) %over deviants
    fname = [fg(mm).folder '/' fg(mm).name];
    %getting MNI coordinates of significant voxels within the provided image
    [ mni_coords, xform ] = osl_mnimask2mnicoords(fname);
    %loading the image
    V = nii.load(fname);
    %extracting statistics
    VV = V(V~=0);
    %indices of non-zero values of nifti image
    VI = find(V~=0);
    %path to AAL template
    parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz'; %load this from the provided codes folder
    %loading AAL labels
    load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat'); %load this from the provided codes folder
    %extracting AAL coordinates information
    K = nii.load(parcelfile);
    %sorting results in order to have strongest voxels at the top (positive t-values) or at the bottom (negative t-values)
    [VV2, II] = sort(VV,'descend');
    VI = VI(II);
    mni_coords = mni_coords(II,:);
    %final cell
    PD = cell(length(VV2),4);
    %getting AAL indices
    ROI = zeros(length(VI),1);
    cnt = 0;
    for ii = 1:length(VI)
        ROI(ii) = K(VI(ii));
        if ROI(ii) > 0 && ROI(ii) < 91
            cnt = cnt + 1;
            PD(cnt,1) = {lab(ROI(ii),3:end)}; %storing ROI
            PD(cnt,4) = {mni_coords(ii,:)}; %storing MNI coordinates
            if mni_coords(ii,1) > 0 %storing hemisphere
                PD(cnt,2) = {'R'};
            else
                PD(cnt,2) = {'L'};
            end
            PD(cnt,3) = {round(VV2(ii),2)}; %storing t-statistics
        end
    end
    PDn = cell2table(PD(~any(cellfun('isempty',PD),2),:)); %remove the possible empty cell
    writetable(PDn,[out '/' fg(mm).name(1:3) '_generalMMN_source.xlsx'],'Sheet',1) %printing excel file
end

%%
