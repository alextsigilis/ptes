import os
import re
import jax
import pickle
import pyedflib
import pandas as pd
import jax.numpy as jnp
import numpy as np
from scipy.ndimage import zoom
from dataloader import load_patient
from hosa import *

# ===================================================================================== #
#                                                                                       #
#                               Script Parameters                                       #
#                                                                                       #
# ===================================================================================== #

# Path to edf reocrdings
datadir = "/data/thmmy/Year_4/Semester_8/PTES/Project/dataset/recordings"

# where to save the estimations
savedir = "./cache/bic"

# How many Segments to use for the estimations
nseg = 32

# The sampling frequncy
fs = 256

# The duration of each epoch
dt = 30

# Frequencies to keep
kl, ku = 150, 250

# ===================================================================================== #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ===================================================================================== #

# Regex for finding the patient files
pattern = re.compile("SN[0-9]{3}\.edf")

# Get all patients in `datadir`
patients = [file[:5] for file in os.listdir(datadir) if pattern.match(file)]


# Create save dir
os.mkdir(savedir)

# Estimating the Bicoherence for every patient
for patient in patients:

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
    estimations = {ch: [] for ch in channels}
    estimations['stage'] = stages

    # For every epoch estimate Bicohernce
    for i,x in enumerate(epochs):
        print(f"\t Epoch {i+1}/{nepochs}", end="\r")
        b = np.array(bicoher(x))
        _, M, M = b.shape
        for j, ch in enumerate(channels):
            estimations[ch].append(b[j,:,:])

    print("\nDone!")
    print(f"Saving to file {patient}.pickle...", end="")

    df = pd.DataFrame(estimations)
    with open(os.path.join(savedir, f"{patient}.pickle"), "wb") as file:
        pickle.dump(df, file)
    print("Done!", end="\n\n")
