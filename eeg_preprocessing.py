"""
========================================================================
 Single-subject ERP pipeline — from raw BDF to single-trial ERPs 
 Author: Ryan Henderson, with help from Claude
========================================================================

Memory note: depending on the size of your EEG file (I'm running a BioSemi 
system so this is a BDF for me), reading the BDF can be taxing on your system
memory (aka RAM). Run this script either on a computer with 16+ GB, or on a cluster.
After unpacking the raw data, we resample immediately, shrinking everything.

"""

import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")      # Use this setting to save PDFs on a headless system
import mne

# ======================================================================
#  PARAMETERS
# ======================================================================

# ---- Files & folders - set these -------------------------------------------------
BDF_PATH   = ""
OUT_ROOT   = ""
SUBJECT_ID = ""

# ---- Trigger / events ------------------------------------------------
# NOTE: Make sure to adapt these to your system's settings if you have a BioSemi
STIM_CHANNEL = "Status"
EVENT_MASK   = 0xFF        # strips BioSemi's 31744 offset (31845 -> 101, etc.)
EVENT_ID     = {f"step{n}": 100 + n for n in range(1, 10)}   # step1->101 … step9->109

# ---- Downsampling ----------------------------------------------------
NEW_SFREQ = 512            # Hz. The original rate is read from the file. We
                           # will use 512 Hz since that is overkill for 30 Hz ERPs. 

# ---- Stimulus-onset timing correction (READ THIS) --------------------
# The trigger marks when the audio FILE started. The sound arrives 51 ms later:
# 50 ms of leading silence in the file + ~1 ms down the Etymotic tubes. To put
# epoch t=0 on TRUE sound onset, each marker moves 51 ms LATER in the recording.
# Zero now sits where the sound began, not where the file started.
#
# SIGN: POSITIVE = sound arrives AFTER the trigger in most cases. Applied in the
# 512 Hz timebase, so 51 ms -> round(0.051 * 512) = 26 samples = 50.78 ms. Redo this
# calculation if you use a different trigger box setup/different onset latency.

STIM_DELAY_MS = 51.0

# ---- Plot dials ------------------------------------------------------
STEPS_TO_PLOT = [1, 2, 3, 4, 5, 6, 7]   # which VOT steps (1–9) go in the ERP figure
PLOT_COMBINE  = "mean"            # average F3/Fz/F4 into one trace per step;
                                  # set to None for one panel per channel

# ---- Channel removal -------------------------------------------------
CHANS_TO_DROP = ["EXG3", "EXG4", "EXG5", "EXG6", "EXG7", "EXG8"]

# ---- Reference -------------------------------------------------------
REF_CHANNELS = ["EXG1", "EXG2"]   # your channels 33 & 34 (linked mastoids)

# ---- Filter ----------------------------------------------------------
L_FREQ = 0.1
H_FREQ = 30.0

# ---- Montage ---------------------------------------------------------
MONTAGE = "biosemi32"

# ---- Artifact rejection (peak-to-peak, per channel type) -------------
# If you recorded from external canthus EEG sensors, use them instead
VEOG_CHANNELS = ["Fp1", "Fp2"]    # vertical EOG proxy (blinks)
HEOG_CHANNELS = ["F7", "F8"]      # horizontal EOG proxy (saccades)
EOG_CHANNELS  = VEOG_CHANNELS + HEOG_CHANNELS
EOG_REJECT_UV = 75                # µV -> eye-monitor channels
EEG_REJECT_UV = 150               # µV -> every other scalp channel

# ---- Epoching --------------------------------------------------------
TMIN     = -0.200
TMAX     =  0.300
BASELINE = (-0.198, 0)            # Set the baseline for displaying ERP 

# ---- Amplitude measurement (window is relative to corrected onset) ---
MEAS_CHANNELS = ["F3", "Fz", "F4"] # Alter to your electrodes of interest
# Set this amplitude to capture the ERP component of your choice
MEAS_TMIN = 0.150                 # s (150 ms post sound onset) 
MEAS_TMAX = 0.180                 # s (180 ms post sound onset)
CSV_PATH  = os.path.join(OUT_ROOT, f"{SUBJECT_ID}_amplitudes.csv")

# Output rescaled events file for transparency/to save rescaling computation
# for future script runs
EVENTS_PATH = os.path.join(OUT_ROOT, "downsampled", f"{SUBJECT_ID}-eve.fif")


# ======================================================================
#  Save file helper
# ======================================================================
def step_path(folder, filename):
    """Create OUT_ROOT/folder/ (like `mkdir -p`) and return the path to filename."""
    d = os.path.join(OUT_ROOT, folder)
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, filename)


# ======================================================================
#  STEP 1 — Read BDF, find events at original rate, resample data+events together
# ======================================================================
raw = mne.io.read_raw_bdf(BDF_PATH, preload=True)
print(f"BDF sample rate read from file: {raw.info['sfreq']:.0f} Hz")

# Extract triggers, stripping the high bits so codes come back as 101..109 (and 
# 14/15, which EVENT_ID ignores).
events = mne.find_events(raw, stim_channel=STIM_CHANNEL,
                         mask=EVENT_MASK, mask_type="and")
print(f"Found {len(events)} events at full rate.")

# Resample data and events simultaneously. MNE reads the true source rate from the 
# object and rescales the event samples to match NEW_SFREQ.
raw, events = raw.resample(NEW_SFREQ, events=events)

# Save the downsampled raw and the canonical (unshifted) events, now correctly
# aligned. A future resume script can load both without any rescaling.
raw.save(step_path("downsampled", f"{SUBJECT_ID}_downsampled_raw.fif"),
         overwrite=True)
mne.write_events(EVENTS_PATH, events, overwrite=True)   # UNSHIFTED triggers, 512 Hz timebase

#  Timing correction — move markers onto true sound onset (in 512 Hz clock)
shift_samples = int(round(STIM_DELAY_MS / 1000.0 * raw.info["sfreq"]))
events[:, 0] += shift_samples           # + = later (toward true sound onset)
print(f"Shifted events by {shift_samples} samples "
      f"({shift_samples / raw.info['sfreq'] * 1000:.2f} ms).")

# Sanity check: how far does the last epoch reach vs. the end of the data?
need_past_end = (events[:, 0].max() + int(TMAX * raw.info["sfreq"])) - raw.last_samp
print(f"Samples the last epoch needs beyond the data end: {need_past_end} "
      f"(<=0 means everything fits; a small positive number = only the final "
      f"trial or two clipped, which is normal).")


# ======================================================================
#  STEP 2 — Remove unused channels
# ======================================================================
if CHANS_TO_DROP:
    raw.drop_channels(CHANS_TO_DROP)
raw.save(step_path("channels_removed", f"{SUBJECT_ID}_chdrop_raw.fif"),
         overwrite=True)


# ======================================================================
#  STEP 3 — Re-reference to channels 33 & 34 (L/R mastoids - change if 
#  you use a different reference)
# ======================================================================
raw.set_eeg_reference(ref_channels=REF_CHANNELS)
raw.set_channel_types({ch: "misc" for ch in REF_CHANNELS})   # keep out of EEG set
raw.save(step_path("rereferenced", f"{SUBJECT_ID}_reref_raw.fif"),
         overwrite=True)


# ======================================================================
#  STEP 4 — Band-pass filter 0.1–30 Hz (continuous data)
# ======================================================================
raw.filter(l_freq=L_FREQ, h_freq=H_FREQ)
raw.save(step_path("filtered", f"{SUBJECT_ID}_filtered_raw.fif"),
         overwrite=True)


# ======================================================================
#  STEP 5 — Add channel locations (montage)
# ======================================================================
montage = mne.channels.make_standard_montage(MONTAGE)
raw.set_montage(montage, on_missing="ignore")


# ======================================================================
#  STEPS 6–8 — Epoch and peak-to-peak artifact rejection 
# ======================================================================
raw.set_channel_types({ch: "eog" for ch in EOG_CHANNELS})

reject = {"eog": EOG_REJECT_UV * 1e-6,     # µV -> V
          "eeg": EEG_REJECT_UV * 1e-6}

epochs = mne.Epochs(raw, events, event_id=EVENT_ID,
                    tmin=TMIN, tmax=TMAX, baseline=BASELINE,
                    reject=reject, preload=True)

epochs.save(step_path("epoched", f"{SUBJECT_ID}_epo.fif"), overwrite=True)

# Report rejections as retained trials vs found trials
print(f"Epochs retained: {len(epochs)} / {len(events)} "
      f"({len(events) - len(epochs)} dropped).")


# =================================================================================
#  STEP 9 — Mean amplitude at F3/Fz/F4, 150–180 ms, single-trial -> long-format CSV
# =================================================================================
rows = []
for cond_name in EVENT_ID:
    ep = epochs[cond_name]
    if len(ep) == 0:
        print(f"WARNING: {cond_name} has 0 surviving trials.")
        continue
    cropped = ep.copy().crop(tmin=MEAS_TMIN, tmax=MEAS_TMAX)
    data = cropped.get_data(picks=MEAS_CHANNELS) * 1e6      # V -> µV
    mean_amp = data.mean(axis=2)                            # mean over time
    step_num = EVENT_ID[cond_name] - 100
    for trial_i in range(mean_amp.shape[0]):
        for ch_i, ch_name in enumerate(MEAS_CHANNELS):
            rows.append({
                "Subject":       SUBJECT_ID,
                "Step":          step_num,
                "Channel":       ch_name,
                "Trial":         trial_i + 1,
                "MeanAmplitude": mean_amp[trial_i, ch_i],
                "WindowMs":      f"{int(MEAS_TMIN*1000)}-{int(MEAS_TMAX*1000)}", # convert from secs -> ms
            })

amp_df = pd.DataFrame(rows)
amp_df.to_csv(CSV_PATH, index=False)

counts = (amp_df[amp_df["Channel"] == MEAS_CHANNELS[0]]
          .groupby("Step").size().reset_index(name="TrialsRetained"))
counts.to_csv(os.path.join(OUT_ROOT, f"{SUBJECT_ID}_trials_per_step.csv"),
              index=False)


# ======================================================================
#  STEP 10 — Average ALL steps into ERPs; PLOT only STEPS_TO_PLOT
# ======================================================================
evokeds = {}
for cond in EVENT_ID:
    if len(epochs[cond]) > 0:
        evokeds[cond] = epochs[cond].average()

mne.write_evokeds(step_path("erps", f"{SUBJECT_ID}_ave.fif"),
                  list(evokeds.values()), overwrite=True)

requested = [f"step{n}" for n in STEPS_TO_PLOT]
plot_evokeds = {name: evokeds[name] for name in requested if name in evokeds}
dropped = [name for name in requested if name not in evokeds]
if dropped:
    print(f"NOTE: requested steps with no data, skipped in the plot: {dropped}")

fig = mne.viz.plot_compare_evokeds(plot_evokeds, picks=MEAS_CHANNELS,
                                   combine=PLOT_COMBINE, show=False)[0]
fig.savefig(step_path("erps", f"{SUBJECT_ID}_erps.pdf"))

print("Done. Steps saved under:", OUT_ROOT)