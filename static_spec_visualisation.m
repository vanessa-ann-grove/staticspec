%% VIPA rsEEG Static Spectral Analysis: Visualisation Code
%
% Vanessa Grove, MSc
%
% The following script is a continuation of the VIPA rsEEG analysis code
% that is specific to the variable 'fft_data'. This code is designed to
% generate several figures as a visualisation tool for scalp-level
% analysis.
%
% To use code: Ensure correct figure options are commented/uncommented based on
% desired analysis. Then update first three lines of code according to data.
%
%% Data Parameters: UPDATE FOR EACH DATASET
clear all

ppt_id = 'V6'; % Participant ID
Exp = 3; % Experiment Number
comparison = 1; %Data Comparison (1 = order, 2 = condition) 


%% Data Setup

% Load Data
data_path = ['S:\VIPA Study\EEG Data\MATLAB\E' num2str(Exp)...
    '_fft_data\\fft_data_' ppt_id '.mat'] %will print to ensure correct dataset has been loaded
load(data_path)

if Exp == 2
    n = 3;
    idn = 1;
else if Exp == 3
    n = 5;
    idn = 2;
    end
end

if comparison == 1
    load('S:\VIPA Study\EEG Data\MATLAB\\experiment_ids_order.mat')
else if comparison == 2
    load('S:\VIPA Study\EEG Data\MATLAB\\experiment_ids_cond.mat')
    end
end

stat1 = n+2; %First row of fft_data statistics
freq_bands = {'Delta' 'Theta' 'Alpha' 'Beta 1' 'Beta 2' 'Beta 3' 'Gamma'};
N = fft_data(1).Npermutes; %number of permutations
chaninfo = fft_data(1).chaninfo;


%% Visualisation

%Figure: Localisation of significant channels on scalp over time
%
%
%
figure(3),clf
figtitle = [ppt_id ': ' 'Significance by channel over time'];

% Define colour limits (opt.)
%
min_max = zeros(2,7);
for freqi = 1:7
    freq_xtrems = zeros(2,n); 
    for datai = stat1:length(fft_data)   
        %dat = fft_data(datai).truestat(freqi,:); %tstat
        dat = fft_data(datai - n).powavg(:,freqi)-fft_data(1).powavg(:,freqi); 
        freq_xtrems(1,(datai-n-1)) = min(dat);
        freq_xtrems(2,(datai-n-1)) = max(dat);
    end
    
    min_max(1,freqi) = min(freq_xtrems(1,:));
    min_max(2,freqi) = max(freq_xtrems(2,:));
end
%}

% Plot
for datai = stat1:length(fft_data) 
    sigData = fft_data(datai).sigdif';
    sigData2 = fft_data(datai).scalpsig';

     % tstat as heatmap data
        %data2plot = fft_data(datai).truestat;
       
     % channel power difference (post-pre) as heatmap data   
        preData = fft_data(1).powavg';
        postData = fft_data(datai - n).powavg';
        data2plot = postData - preData;
    
    for freqi = 1:size(data2plot,1)
        mi = min_max(1,freqi);
        ma = min_max(2,freqi);
        clear sigchanlocs
        for chani = 1:64
            if sigData(chani,freqi)== 1
                sigchanlocs(chani) = chaninfo(chani);
            end
        end
   %
        if exist('sigchanlocs', 'var') %&& sigData2(freqi)==1
            subplot(n,7,(freqi+(7*(datai-(n+2)))))
            topoplot(data2plot(freqi,:),chaninfo,'electrodes','off','style','map','maplimits',[mi ma]);
            hold on
            topoplot([],sigchanlocs,  'electrodes', 'on');
            colorbar
            hold off
        else
        %}
            subplot(n,7,(freqi+(7*(datai-(n+2)))))
            topoplot(data2plot(freqi,:),chaninfo, 'electrodes','off','style','map','maplimits',[mi ma]);
            colorbar
        end
    end
end
set(gcf, 'name', figtitle)
%}

%Figure: Bar graph of whole scalp data with full spectrum 
%
figure(5),clf
    %Compare two sets of data
    %{
combined_dat5 = [fft_data(d1).scalpavg;fft_data(d2).scalpavg]'; 
bar(combined_dat5)
hold on
for freqi=1:size(data2plot2,2)
        if sigData2(freqi) == 1 %>=thresh
           plot([(freqi-0.4),freqi+0.4], [1 1] + max(combined_dat5(freqi,1),combined_dat5(freqi,2)), '-k', 'LineWidth',2)
           plot(mean([(freqi-0.4),freqi+0.4]), 1.3+max(combined_dat5(freqi,1),combined_dat5(freqi,2) ),'*k')
        end
end
%}
    %Compare all data
subplot(212)
vals = experiment_ids(idn).vals;
combined_dat5 = [fft_data(1).scalpavg;fft_data(2).scalpavg;fft_data(3).scalpavg;fft_data(4).scalpavg...
    ;fft_data(5).scalpavg;fft_data(6).scalpavg]';
b = bar(combined_dat5);
b(1).FaceColor = 'k';
hold on
for datai = stat1:length(fft_data)
    sigData2 = fft_data(datai).scalpsig';
    for freqi=1:7
        maxval = combined_dat5(freqi,(datai-n));
        if sigData2(freqi) == 1
           plot(freqi+vals(datai-n), 1.1*maxval ,'*k')  
        end
    end
end
xticklabels({'Delta (1-4 Hz)' 'Theta (4-8 Hz)' 'Alpha (8-12 Hz)'...
    'Beta 1 (12-16 Hz)' 'Beta 2 (16-23 Hz)' 'Beta 3 (23-30 Hz)' 'Gamma (30-45 Hz)'})
ylabel('Average Power (\muV^2)')
if comparison == 1 && Exp == 3
    legend({'Baseline' 'G1' 'G2' 'G3' 'G4' 'Post-VR' 'Significant to baseline'})
else if comparison == 2 && Exp == 3
    legend({'Baseline' 'AS' 'ES' 'FF' 'PP' 'Post-VR' 'Significant to baseline'})
else if comparison == 1 && Exp == 2
    legend({'Baseline' 'E1' 'E2' 'Post-VR' 'Significant to baseline'})
        end
    end
end

title([ppt_id ': ' 'Average Scalp Power'])

% Full Spectra
hold on
freq_bounds = dsearchn(fft_data(1).frex, [2 45]');
frex = fft_data(1).frex(freq_bounds(1):freq_bounds(2));
subplot(211)
for i = 1:length(experiment_ids(idn).data)
    powvals = fft_data(i).allpower(:,freq_bounds(1):freq_bounds(2)); 
    spectra = mean(powvals,1);
    p = plot(frex,spectra);
    if i == 1
       p.Color = 'k';
       p.LineWidth = 1.5;
    end
    hold on
end
ax = gca;
ax.XLim = [2 45];
ax.XTick = [2 4 8 12 16 23 30 45];
for o = 1:length(ax.XTick)
    xline(ax.XTick(o))
end
if comparison == 1 && Exp == 3
    legend({'Baseline' 'G1' 'G2' 'G3' 'G4' 'Post-VR'})
else if comparison == 2 && Exp == 3
    legend({'Baseline' 'AS' 'ES' 'FF' 'PP' 'Post-VR'})
else if comparison == 1 && Exp == 2
    legend({'Baseline' 'E1' 'E2' 'Post-VR'})
        end
    end
end
title('Spectra of Controls During VR Intervention')
xlabel('Frequency (Hz)')
ylabel('Log Power Spectral Density 10*log_{10}(\muV^{2}/Hz)')

%}

%Figure: Topo plots- change in power over time for each frequency band
%
% Define map colour limits
min_max = zeros(2,7);
for freqi = 1:7
    freq_xtrems = zeros(2,n+1);
    for datai = 1:(n+1)
        pow = fft_data(datai).powavg(:,freqi);
        freq_xtrems(1,datai) = min(pow);
        freq_xtrems(2,datai) = max(pow);
    end
    min_max(1,freqi) = min(freq_xtrems(1,:));
    min_max(2,freqi) = max(freq_xtrems(2,:));
end

% Create figure with maps
figtitle = [ppt_id ': ' 'Average Scalp Power Over Time'];
figure(20),clf
hold on
for freqi = 1:7
    mi = min_max(1,freqi);
    ma = min_max(2,freqi);
    for datai = 1:(n+1)       
        pow = fft_data(datai).powavg(:,freqi);
        
        subplot((n+1),7,(freqi+(7*(datai-1))))
        topoplot(pow,chaninfo,'electrodes','off','maplimits',[mi ma]);
        colorbar
    end
end
set(gcf, 'name', figtitle)
%}

% Figure: Full Spectral Image (Standalone)
%{
freq_bounds = dsearchn(fft_data(1).frex, [2 45]');
frex = fft_data(1).frex(freq_bounds(1):freq_bounds(2));
figure(2), clf

for i = 1:6
    powvals = fft_data(i).allpower(:,freq_bounds(1):freq_bounds(2)); 
    spectra = mean(powvals,1);
    p = plot(frex,spectra);
    if i == 1
       p.Color = 'k'
       p.LineWidth = 1.5;
    end
    hold on
end
ax = gca;
ax.XLim = [2 45];
legend({'Baseline' 'G1' 'G2' 'G3' 'G4' 'Post-VR'});
title('Spectra of Pain Patients During VR Intervention')
xlabel('Frequency (Hz)')
ylabel('Log Power Spectral Density 10*log_{10}(\muV^{2}/Hz)')
%}

%% Code for older figures
%{
%Figures: Topographical plots for each significant frequency band
%{
for freqi = 1:7
    if sigData2(freqi) == 1    
        figtitle = freq_bands{freqi};
        pre_pow = fft_data(d1).powavg(:,freqi);
        post_pow = fft_data(d2).powavg(:,freqi);
        stat = fft_data(s3).truestat(freqi,:);
        dif = post_pow - pre_pow;
    
        figure(freqi+20),clf
        hold on
        subplot(221)
            topoplotIndie(pre_pow,chaninfo);
            title('Pre-Intervention Data')
            colorbar
        subplot(222)
            topoplotIndie(post_pow,chaninfo);
            plot_title = fft_data(d2).id;
            title(plot_title)
            colorbar
        subplot(223)
            topoplotIndie(dif,chaninfo);
            title('Difference (post-pre)')
            colorbar
        subplot(224)
            topoplotIndie(stat,chaninfo);
            title('Paired t statistic')
            colorbar
        set(gcf, 'name', figtitle)
    end
end
%}

%Figure: Visualising power spectra and average values from single channel
%{
comp2plot = 15;
fig2data1 = fft_data(1).allpower;
fig2data2 = fft_data(2).allpower;

figure(2),clf
figtitle = ['Channel' num2str(comp2plot)];
subplot(311)
plot(hz,fig2data1(comp2plot,1:length(hz)),'k-','linewidth',1.5)
xlabel('Frequency (Hz)')
ylabel('Power (\muV)')
title('FFT Power spectrum: Pre-Intervention')
set(gca,'xlim',[2 40])
subplot(312)
plot(hz,fig2data2(comp2plot,1:length(hz)),'k-','linewidth',1.5)
xlabel('Frequency (Hz)')
ylabel('Power (\muV)')
title('FFT Power spectrum: Post-Intervention')
set(gca,'xlim',[2 40])
subplot(313)
combined_dat = [fft_data(1).powavg(comp2plot,:);fft_data(2).powavg(comp2plot,:)]';
bar(combined_dat)
xticklabels({'delta' 'theta' 'alpha' 'beta 1' 'beta 2' 'beta 3' 'gamma'})
ylabel('Average Power (\muV)')
legend('Pre Intervention', 'Post Intervention')
title(['Average Power Channel' num2str(comp2plot)])
set(gcf, 'name', figtitle)
%}

%Figure: Heatmap with outlined significance
%{
figure(3),clf
subplot(3,1,[1,2])
imagesc(data2plot)
axis xy
hold on
xtk = get(gca, 'XTick');                               
set(gca, 'XTick', [1:7], 'XTickLabel', freq_bands)
for cmpti=1:size(data2plot,1)
    for freqi=1:size(data2plot,2)
        if sigData(cmpti,freqi) == 1 %>=thresh
           rectangle('Position', [(freqi-0.5) (cmpti-0.5) 1 1], 'LineWidth',1.5)
        end
    end
end
axis square
ylabel('Channel Number')
title({['Paired ttest stat with outlined significance'] ['N =' num2str(N) 'permutes']}) 
colorbar

subplot(313)
imagesc(data2plot2)
axis xy
hold on
xtk = get(gca, 'XTick');                               
set(gca, 'XTick', [1:7], 'XTickLabel', freq_bands)
    for freqi=1:size(data2plot2,2)
        if sigData2(freqi) == 1 %>=thresh
           rectangle('Position', [(freqi-0.5) 0.5 1 1], 'LineWidth',1.5)
        end
    end
xlabel('Frequency Band') 
title('Scalp Averaged Data')
colorbar
%}
%}
