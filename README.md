The present script performs automated EEG preprocessing for single-trial ERP-based analyses. The dataset for which this script was designed contained eventcodes already, so there is no insertEventcodes step. If you have a dataset without eventcodes, make sure to insert them before running this script. 
In order, the analysis is the following:
1) Set important variables
    - EVENTID = conversion of BioSemi eventcodes to more human-readable codes
    - STIM-DELAY = 50 ms of silence prior to each stimulus onset, plus ~1 ms to account for the delay in the sound travelling from the earbud speaker to the subject's ears
    - L_freq, H_freq = bandpass filter frequencies that are optimal for auditory ERP research. See https://doi.org/10.1111/psyp.12437Digital 
2) Read in the data, downsample to 512 Hz
    - BioSemi systems tend to deafult to very high sampling rates. This makes very large files that are slow to process. We will downsample to 512 Hz, because that is well above the Nyquist frequency for capturing ERP effects. 
3)  Remove unused channels (the recording montage only used the first two EXG channels, so we'll drop the rest)
4) Re-reference
5) Bandpass filter the EEG signal 
6) Add channel locations 
    - Always a good step in case someone would like to do source localization in the future
7) Epoch data
8) Artifact Rejection
9) Output average EEG amplitude at the N100 latency
10) Create averaged ERPs, and output a test figure for quality control
