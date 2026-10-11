# eeg_preprocessing

A single-subject EEG preprocessing pipeline that takes a continuous BioSemi recording from raw BDF to averaged ERPs and single-trial amplitude measurements.

![Averaged ERPs across VOT steps, measured at F3/Fz/F4](example_erps.png)

Your figure will not look identical — `mne.viz` styling varies between versions — but the general shape should resemble this if you are looking at N100 data.

## Requirements

    pip install -r requirements.txt

Tested on Linux-4.18.0, Python 3.13.5, MNE 1.10.0, NumPy 2.2.6, pandas 3.0.5, Matplotlib 3.10.3 (backend=Agg).

Reading a full-rate BioSemi BDF is memory-hungry. Run this on a machine with 16+ GB of RAM or on a cluster. The script resamples immediately after loading, so everything downstream is much smaller.

## Usage

Set the variables in the `PARAMETERS` block at the top of `eeg_preprocessing.py`, then:

    python eeg_preprocessing.py

| Variable | Set it to |
| --- | --- |
| `BDF_PATH` | Path to the raw BioSemi `.bdf` file |
| `OUT_ROOT` | Directory where all outputs are written |
| `SUBJECT_ID` | Identifier used to name every output file |
| `EVENT_ID` | Mapping from your trigger codes to condition names |
| `STIM_DELAY_MS` | Delay between the trigger and true sound onset at the ear |
| `CHANS_TO_DROP` | Unused channels to discard |
| `REF_CHANNELS` | Offline reference channels (linked mastoids by default) |
| `MONTAGE` | MNE montage name matching your cap |
| `VEOG_CHANNELS` / `HEOG_CHANNELS` | Channels used to detect blinks and saccades |
| `MEAS_CHANNELS` / `MEAS_TMIN` / `MEAS_TMAX` | Electrodes and time window for amplitude measurement |

This script assumes your recording already contains event codes. If it does not, insert them before running.

## Inputs and outputs

**Input:** one continuous BioSemi `.bdf` recording, with event codes on the `Status` channel.

**Outputs**, all written under `OUT_ROOT`:

| Path | Contents |
| --- | --- |
| `downsampled/` | Resampled raw data, plus the unshifted event file on the new timebase |
| `channels_removed/`, `rereferenced/`, `filtered/` | Intermediate `.fif` files after each stage |
| `epoched/` | Epochs surviving artifact rejection |
| `erps/*_ave.fif` | Averaged ERP for each condition |
| `erps/*_erps.pdf` | Quality-control figure |
| `*_amplitudes.csv` | Single-trial mean amplitudes, long format |
| `*_trials_per_step.csv` | Trials retained per condition after rejection |

Saving at every stage means a failed run can be picked up partway through, and makes it possible to find the stage where a problem was introduced.

## How it works

1) Set important variables
    - EVENTID = conversion of BioSemi eventcodes to more human-readable codes.
    - STIM-DELAY = 50 ms of silence prior to each stimulus onset, plus ~1 ms to account for the delay in the sound travelling from the Etymotic ER-3A earbud speakers through tubes to the subject's ears.
    - L_freq, H_freq = bandpass filter frequencies that are optimal for auditory ERP research. See https://doi.org/10.1111/psyp.12437 for optimal filtering suggestions.
2) Read in the data, downsample to 512 Hz.
    - BioSemi systems record at very high sampling rates. This makes very large files that are slow to process. We will downsample to 512 Hz, because that is well above the Nyquist frequency for capturing 30 Hz ERPs.
3) Remove unused channels (the recording montage only used the first two EXG channels, so we'll drop the rest).
4) Re-reference the continuous signal to your choice of offline reference (I prefer linked mastoids).
5) Bandpass filter the EEG signal, removing slow sweat artifacts along with high-frequency artifacts.
6) Add channel locations
    - Always a good step in case someone would like to do source localization in the future
7) Epoch data, converting the long continuous signal into chunks that begin 200 ms before stimulus onset and end 300 ms post-onset. Adjust the second time bound if you're investigating ERPs at a different latency.
8) Artifact Rejection - reject trials that are contaminated by blinks and/or muscle movements.
9) Output average EEG amplitude on each trial over the chosen latency.
10) Create averaged ERPs, collapsing across trials, and output a test figure for quality control.

## License

MIT — see [LICENSE](LICENSE).
