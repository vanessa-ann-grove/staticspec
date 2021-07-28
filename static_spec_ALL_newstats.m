%% VIPA rsEEG Static Spectral Analysis- ALL TIMEPOINTS
%
% Vanessa Grove, MSc
% 
% The following script is designed for the basic analysis of resting state
% electrophysiogical data. For use on EEGLab-formatted, preprocessed data
% which has been properly inspected and cleaned of artefacts.
%
% To use code: Update first three lines of code according to desired analysis. 
% If new analysis, keep code as is. If repeat analysis for existing data,
% comment out sections 3 and 4 and uncomment first part of section 5. 
%
% Note: Only will work on study laptop, when connected to UEA server.

%% Section 1: Data Parameters- UPDATE FOR EACH DATASET!

ppt_id = 'R10'; %Participant ID
Exp = 3; %Experiment Number (2 or 3)
comparison = 1; %Data Comparison (1 = order, 2 = condition (group stats only)) 
rerun_statistics = 'Y'; %{'Y' if rerunning statistics for existing fft_data variable,...
                        % 'N' if running EEG data to create new
                        % fft_data variable}

% If rerunning statistics, adjust parameters here
test_stat = 'Paired Wilcoxon Signed Rank';
uncorrected_pval = .01;
critical_pval = uncorrected_pval / 7; %Bonferonni correction for 7 freq bands
N = 5000; %number of permutes


%% Section 2: Data Setup

if Exp == 2
    n = 3;
    idn = 1;
    filepath = 'S:\\VIPA Study\\EEG Data\\Experiment 2\\Resting State Data\\Processed Files\\';
    if comparison == 1
        condition = {'pre_EC' 'E1' 'E2' 'post_EC'};
    else if comparison == 2
         condition = {'pre_EC' 'sun' 'sno' 'post_EC'};
        end
    end
else if Exp == 3
    n = 5;
    idn = 2;
    filepath = 'S:\\VIPA Study\\EEG Data\\Experiment 3\\Resting State Data\\Processed Files\\';
    if comparison == 1
         condition = {'EC_pre' 'G1' 'G2' 'G3' 'G4' 'EC_post'};
    else if comparison == 2
         condition = {'EC_pre' 'AS' 'ES' 'FF' 'PP' 'EC_post'};
        end
    end
end
end

if comparison == 1
    load('S:\VIPA Study\EEG Data\MATLAB\\experiment_ids_order.mat')
else if comparison == 2
    load('S:\VIPA Study\EEG Data\MATLAB\\experiment_ids_cond.mat')
    end
end
data_ids = experiment_ids(idn).data;

%% Section 3: Data Load

if rerun_statistics == 'N'
%Load preprocessed EEG files into an ALLEEG variable structure

all_files = cell.empty(length(condition),0);
for time = 1:length(condition)
    file = [char(ppt_id), '_' char(condition(time)) '.set'];
    all_files{time} = char(file);
end

clear EEG ALLEEG CURRENTSET ALLCOM
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename',all_files,'filepath', filepath);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'study',0);


%% Section 4: EEGLab-based Power Calculations

% Initialize Fourier output matrix
fft_data = struct([]);
fft_data(1).Ppt_ID = ppt_id;
fft_data(1).Exp = Exp;
fft_data(1).Test_Stat = test_stat;
fft_data(1).Critical_Pval = critical_pval;
fft_data(1).Npermutes = N;
fft_datalabels = experiment_ids(idn).fft;

for j=1:length(ALLEEG)
  % Specify Data
    data_struct = ALLEEG(j);  
    data = data_struct.data; 
    nData = length(data_struct.times); %number of data points
    srate = data_struct.srate; %sampling rate
    nbchan = size(data,1); %number of channels
    
  % Store Data
    fft_data(j).id = fft_datalabels{j}; 
    fft_data(j).chaninfo = ALLEEG(j).chanlocs; %channel information
    fft_data(j).voltagedata = data; %store voltage data in output variable
    
    %Write voltage data to .csv file for use in sLORETA (opt.)
    %{
    id = data_ids{j};
    name = [id '_' ppt_id '.csv'];
    path = ['C:\\Users\\The UEA VR & EEG Lab\\Documents\\sLORETA\\EEG Data\\Experiment ' num2str(Exp) '\\' name]; 
    pop_export(data_struct,path,'elec','off','time','off','precision',4);
    %}
    
    %Compute FFT for each channel (will return power vals in dB)
    for i=1:nbchan
    figure(2),clf
    [spectra,freqs] = spectopo(data(i,:,:), 0, srate);
    fft_data(j).frex = freqs;
    fft_data(j).allpower(i,:) = spectra;

        %Freq Band Indices
        deltaIdx = find(freqs>=2 & freqs<=4);
        thetaIdx = find(freqs>=4 & freqs<=8);
        alphaIdx = find(freqs>=8 & freqs<=12);
        beta1Idx  = find(freqs>=12 & freqs<=16);
        beta2Idx  = find(freqs>=16 & freqs<=23);
        beta3Idx  = find(freqs>=23 & freqs<=30);
        gammaIdx = find(freqs>=30 & freqs<=45);

        % compute absolute power (unit: uV^2/Hz)
        fft_data(j).powavg(i,1) = 10^(mean(spectra(deltaIdx))/10);
        fft_data(j).powavg(i,2) = 10^(mean(spectra(thetaIdx))/10);
        fft_data(j).powavg(i,3) = 10^(mean(spectra(alphaIdx))/10);
        fft_data(j).powavg(i,4) = 10^(mean(spectra(beta1Idx))/10);
        fft_data(j).powavg(i,5) = 10^(mean(spectra(beta2Idx))/10);
        fft_data(j).powavg(i,6) = 10^(mean(spectra(beta3Idx))/10);
        fft_data(j).powavg(i,7) = 10^(mean(spectra(gammaIdx))/10);
    end
    fft_data(j).scalpavg = mean(fft_data(j).powavg,1);
end

else if rerun_statistics == 'Y'
% Load existing data (if re-reunning statistical analysis)
data_path = ['S:\VIPA Study\EEG Data\MATLAB\E'...
    num2str(Exp) '_fft_data\\fft_data_' ppt_id '.mat'] 
load(data_path)
%Update info
fft_data(1).Test_Stat = test_stat;
fft_data(1).Critical_Pval = critical_pval;
fft_data(1).Npermutes = N;
    end
end

%% Section 5: Statistics Setup

stats_id = experiment_ids(idn).stats;
d1 = 1; % Baseline dataset (1 = pre-intervention)

%Frequency band boundries
hz = fft_data(d1).frex;
freq_bounds = dsearchn(hz, [2 4 8 12 16 23 30 45]');

for d2 = 2:(length(stats_id)+1)
%Specify output location
s1 = d2+n; 
fft_data(s1).id = stats_id{d2-1}; 

%Data
preInt = fft_data(d1).allpower; 
postInt = fft_data(d2).allpower;
ncmp = min(size(fft_data(d1).powavg,1), size(fft_data(d2).powavg,1)); %number of channels

%figure(d2), clf, hold on
for k=2:length(freq_bounds)
%% Section 6: Statistical Analysis (Whole Scalp: concatenated data)
        %Step 1: Select data from freq band
        preData  = preInt(:,freq_bounds(k-1):freq_bounds(k));
            n1 = size(preData,2);
        postData = postInt(:,freq_bounds(k-1):freq_bounds(k));
            n2 = size(postData,2);
        allData = horzcat(preData,postData)';
        truelabels = cat(1,ones(n1,1),2*ones(n2,1));
    
        % Step 2: Compute true test statistic 
            % Wilcoxon signed rank test
                [p h stats] = signrank(mean(allData(truelabels==1,:),1),...
                    mean(allData(truelabels==2,:),1));
                true_stat = p;
                fft_data(s1).scalpstat(k-1,:) = true_stat;
            
        % Step 3: Permute
        for permi=1:N
            % Step 4: Shuffle data labels
            shuflabels = truelabels(randperm(n1+n2));
            %Step 5: Compute test stat with shuffled labels and store
                [p h stats] = signrank(mean(allData(shuflabels==1,:),1),...
                    mean(allData(shuflabels==2,:),1));
                fft_data(s1).scalpperm(k-1,permi) = p;
        end
        
        %Step 6: Calculate p value and create new variable to only store significant test stats
        scalp_stat = fft_data(s1).scalpstat(k-1);
        abs_xtreme = fft_data(s1).scalpperm(k-1,:);
        pval = sum(abs_xtreme < scalp_stat,2) / N;
        fft_data(s1).scalppvals(k-1) = pval;
        trialsig = pval  < critical_pval;
        fft_data(s1).scalpsig(k-1) = trialsig; 

%Plotting (opt.)
%{
labels = {'Max Perm' 'Delta' 'Theta' 'Alpha' 'Beta 1' 'Beta 2' 'Beta 3' 'Gamma'};
    subplot(7,1,(k-1))
    hold on
    histogram(fft_data(s1).scalpperm(k-1,:))
    plot([1 1]*fft_data(s1).scalpstat(k-1),get(gca,'ylim'),'r--','linew',2)
    xlabel('Mean value')
    ylabel('Count')
    %title(stats_id{d2-1})
    %legend(labels)
    figtitle = [ppt_id, '- Stats'];
    set(gcf, 'name', figtitle)
%}

%% Section 7: Statistical Analysis (Single channel: concatenated data)
% Only to be performed if freq band is significant at scalp level
%
trialsig = fft_data(s1).scalpsig(k-1);
if trialsig == 0
        fft_data(s1).truestat(k-1,:) = 1+zeros(1,ncmp);
        fft_data(s1).chanperm(k-1,:) = zeros(1,N);
        fft_data(s1).chan_pvals(k-1,:) = 1+zeros(1,ncmp);
        fft_data(s1).sigdif(k-1,:) = zeros(1,ncmp);    
else if trialsig == 1
        % Step 7: Compute true test statistic 
            %Wilcoxon Signed Rank test for paired samples
            for chani = 1:ncmp
                [p h stats] = signrank(allData(truelabels==1,chani),allData(truelabels==2,chani));
               fft_data(s1).truestat(k-1,chani) = p;  %
            end
        % Step 8a: Permute and store max value
        for permi=1:N
        % Step 8b: Shuffle data labels
            shuflabels = truelabels(randperm(n1+n2));
        %Step 8c: Compute test stat with shuffled labels and store
                perm_stat = zeros(1,ncmp);
                for chani = 1:ncmp
                    [p h stats] = signrank(allData(shuflabels==1,chani),allData(shuflabels==2,chani));
                    perm_stat(chani) = p;
                end
                 fft_data(s1).chanperm(k-1,permi) = min(perm_stat);
        end
        %Step 9: Calculate p value and create new variable to demonstrate significant results
        chan_stat = fft_data(s1).truestat(k-1,:);
        abs_xtreme =  fft_data(s1).chanperm(k-1,:)';
        pval = sum(abs_xtreme < chan_stat) / N;
        fft_data(s1).chan_pvals(k-1,:) = pval;
        trialsig = pval < critical_pval;
        fft_data(s1).sigdif(k-1,:) = trialsig;
            end
        end
    end
end

%% Print Results for User View

disp('Analysis Complete. Dont forget to save fft_data variable!!!!')

statistics = experiment_ids.data;
freq_bands = {'Delta band' 'Theta band' 'Alpha band' 'Beta1 band'...
    'Beta2 band' 'Beta3 band' 'Gamma band'};

if Exp == 3
all_sig = vertcat(fft_data(7).scalpsig,fft_data(8).scalpsig,...
    fft_data(9).scalpsig,fft_data(10).scalpsig,fft_data(11).scalpsig);
else if Exp == 2
  all_sig = vertcat(fft_data(5).scalpsig,fft_data(6).scalpsig,...
    fft_data(7).scalpsig);  
    end
end

if any(all_sig,'all') == 1
    disp('Significant Results Observed:')
for s = 1:(length(statistics)-1)
    for f = 1:length(freq_bands)
       if all_sig(s,f) == 1 
           disp([char(statistics(s+1)) ' ' char(freq_bands(f))...
               ' is significant to baseline!'])
       end
    end
end
else if any(all_sig, 'all') == 0
        disp('No significant results observed')
    end
end
