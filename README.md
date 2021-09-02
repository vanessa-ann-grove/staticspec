These codes are specific to data collected for the VIPA Study at the University of East Anglia, led by Dr Jordan Tsigarides. Written by Vanessa Grove (MSc)

Purpose of codes: Statistics code runs spectral analysis on resting state EEG data and compare data during/following use of Virtual Reality technology to a baseline dataset. Main analysis occurs using data averaged across the whole scalp in seven frequency bands. If the power in a frequency band is found to be significantly different to the baseline dataset, a supplementary analysis will run which will analyse each of the 64 channels individually to examine which channels/groups of channels are demonstrating differences. The results of this code can be visualised in the visualisation code.

IMPORTANT- These codes can only be run on a device that is connected to the S: drive on a UEA-based server.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STATISTICS CODE (static_spec_ALL_newstats.m)

Setup: To use this code on new data, user will update the first section of code with the appropriate information. This must be done manually before running the code.

Code Structure: 
SECTION 1- Data Parameters

* ppt_id- ID number of participant or appended dataset. IMPORTANT: this variable must match the processed EEG data file of interest exactly (i.e. if EEG file is called APPENDED_EC_pre, then 'ppt_id' MUST be 'APPENDED').
* Exp- Number of experiment being analysed (either 2 or 3)
* Comparison- Type of comparison of interest (1 = sequential order (default), 2 = condition (i.e. by game or by environment).
* rerun_statistics- 'Y' if rerunning statistics with different parameters on an existing fft_data variable (default). 'N' if creating a new fft_data variable from preprocessed EEG files

Other information in this section regarding statistical analysis details will get stored in the output variable for future reference. Update this section if rerunning statistical analysis with new parameters.

SECTION 2- Data Parameters This section will set up the filepath for accessing the EEG files of interest and will use the information provided by the first section to prepare the output variable according to the experiment number and type of comparison being facilitated.

SECTION 3- Data Load This section will load preprocessed EEG files of interest into an 'ALLEEG' structure if creating a new fft_data variable. Otherwise, the existing fft_data variable will be loaded in from memory and section 4 will be skipped.

SECTION 4- EEGLab-based Power Calculations
* Output variable 'fft_data' is initialized and set up
* Each of the EEG datasets loaded into the 'ALLEEG' variable will be transformed from voltage to power data using a Fourier Transform. This is achieved using the EEGLab command ‘spectopo’. The output of this calculation is a 64 x 126 matrix, which represents the absolute power (in decibels, average over all epochs) at each frequency 0-125 Hz for each of the 64 channels. This is stored in the output variable under 'fft_data(n).allpower' where n = 1:6 for experiment 3 or n = 1:4 for experiment 2.
* The data from this matrix is then grouped into the desired frequency bands. Absolute power for each frequency band is calculated by averaging the power value of each of the frequencies within the range of each band. This number is then converted to microvolts and stored in 'fft_data(n).powavg'.
* Finally, the microvolt data is averaged across all 64 channels to create a single scalp value. This is stored in 'fft_data(n).scalpavg'.

SECTION 5- Statistics Setup This section of code prepares the output variable to store the results of the statistical analysis. The ‘for’ loop that begins in this section cycles through a comparison of the baseline dataset to all of the following datasets (i.e. baseline:G1, baseline:G2… baseline:post-intervention).

SECTION 6- Whole-Scalp Analysis
* ‘For’ loop cycles through each of the seven frequency bands of interest
* Data is selected from fft_data(1).allpower and fft_data(n).allpower for each of the frequencies included within the frequency range (output: 64 x m matrices, where m = number of frequencies in specified band at a 1Hz resolution).
* The average across channels is computed for each of the datasets (output: 1 x m matrices)
* The ‘true’ test statistic is calculated via a Wilcoxon Signed Rank test. The p-value of this test is stored as the test statistic and is what is used to build the data-based distribution for permutation testing.
* The data is permuted by shuffling the data points randomly between the two datasets 5000 times and performing a Wilcoxon Signed Rank test for each permute. The output of this step is a 1 x 5000 matrix in which each entry corresponds to a p-value for the shuffled data.
* A final p-value for the permutation testing is calculated by summing the number of times permuted test statistic values were less than the true test statistic, divided by the total number of permutes: (∑(permuted statistic < test statistic)/N). This value is stored in decimal form in 'fft_data.scalppvals'. Additionally, if the p-value is less than the Bonferroni-corrected critical value, this information is stored as a ‘1’ in 'fft_data.scalpsig' and if the value is greater than the critical p-value, it is stored as a ‘0’.
* Optional visualization: The null hypotheses distribution can be visualized as histograms with the true test statistic plotted on top.

SECTION 7- Single Channel Analysis Note: This section is only performed if the result of the whole-scalp data is significant.
* The relevant value from 'fft_data.scalpsig' is evaluated. If the value is ‘0’, this section is skipped.
* The same procedure as above is repeated, however, this time it is computed for each of the 64 channels individually, therefore the ‘true’ test statistic output is a 1 x 64 matrix of p-values corresponding to each channel. The result of permutation testing is a 5000 x 64 matrix. In order to calculate the final p-value for each of the channels and to account for multiple comparisons, the minimum value of each of the permutes is selected (across channels), resulting in a final distribution of 5000 p-values. The true test statistic of each of the 64 channels is then compared to this single distribution using the formula above. The values are stored in 'fft_data.chan_pvals', and the logical values representing significance are stored in 'fft_data.sigdif'.

IMPORTANT- 'fft_data' VARIABLE MUST BE **MANUALLY SAVED** FOLLOWING THE RUNNING OF THIS CODE!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VISUALISATION CODE (static_spec_visualisation.m)

Setup: To use this code on new data, user will update the first four lines of code with the appropriate information. This must be done manually before running the code. Resulting figures (arbitrary numbers) are as follows:

* Figure 20: Average scalp power in each frequency band (columns) at each time point (rows). 
* Figure 5: Top panel displays the full spectra at each time point. Bottom pannel shows the average band power at each time point in bar graph format for easy comparison. Significant differences are displayed with stars over each bar that is different to the baseline measurement.
* Figure 3: Each row of scalp maps represents the mean power difference of each dataset compared to baseline (i.e. row 1 is G1 - baseline for experiment 3 or E1 - baseline for experiment 2) in each of the seven frequency bands (columns). If a significant difference was found at scalp level, channels that were found to be significant will be displayed over the scalp. 
