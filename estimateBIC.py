import os
import re
import jax
import pickle
import pyedflib
import pandas as pd
import jax.numpy as jnp
import numpy as np
from collections import defaultdict
from scipy.ndimage import zoom
from dataloader import load_patient
from hosa import *

# ===================================================================================== #
#                                                                                       #
#                               Script Parameters                                       #
#                                                                                       #
# ===================================================================================== #

# Path to edf reocrdings
datadir = "./physionet.org/files/hmc-sleep-staging/1.0.0/recordings"

# where to save the estimations
savedir = "./cache"

# How many Segments to use for the estimations
nseg = 16

# The sampling frequncy
fs = 256

# The duration of each epoch
dt = 30

# Define for which patients to estimate the bicoherency
start, stop = 0, 70

# ===================================================================================== #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ===================================================================================== #

# Regex for finding the patient files
pattern = re.compile("SN[0-9]{3}\.edf")

# Get all patients in `datadir`
patients = [file[:5] for file in os.listdir(datadir) if pattern.match(file)]
patients = sorted(patients)

# Estimating the Bicoherence for every patient
for pid, patient in enumerate(patients[start:stop]):

    print(f"Loading patient {patient}...", end="")

    # Load data for patient
    signals, headers, _ = pyedflib.highlevel.read_edf(os.path.join(datadir, f"{patient}.edf"))
    scores = pd.read_csv(os.path.join(datadir, f"{patient}_sleepscoring.txt"))
    # Get the names of the channels
    channels = [header['label'] for header in headers]
    print("Done.")
    print("Preprocessing data...", end="")

    # Find all usefull annotations
    idx = scores[' Annotation'].apply(
        lambda s: re.fullmatch("Sleep.*", s.strip()) is not None
    )
    scores = scores[idx]

    # Get only the stages ('W' or 'Ni' or 'R', i=1,2,3)
    stages = [
        s.split(" ")[-1] for s in scores[' Annotation']
    ]

    # Get dimentions of data
    nchan, _ = signals.shape
    nepochs, _, = scores.shape

    # Number of samples per epoch
    N = int(fs*dt)

    # Split into epochs
    epochs = [
        signals[:, i*dt*fs:(i+1)*dt*fs].reshape((nchan, nseg, N//nseg))
        for i in range(nepochs)
    ]

    print("Done.")
    print(f"Estimating Bicoherence for patients {patient}...")

    # Create empty dictionary to store the results
    estimations = defaultdict(list)
    estimations['stage'] = stages

    # For every epoch estimate Bicohernce
    for i,x in enumerate(epochs):
        print(f"\t Epoch {i+1}/{nepochs}", end="\r")
        b = np.array(bicoher(x)).astype(np.float32)
        _, M, M = b.shape
        for j, ch in enumerate(channels):
            estimations[ch].append(b[j,:,:])

    print("\nDone!")

    print(f"Saving to file {patient}.pickle...", end="")
    pd.DataFrame(estimations).to_pickle(f"{savedir}/bic_{patient}.pickle")
    print("Done!")

print("\t Goodbye")

