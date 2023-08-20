import pyedflib
import numpy as np
import pandas as pd
import re

def load_patient(edf, scores, dt=30, fs=256):
    signals, headers, _ = pyedflib.highlevel.read_edf(edf)
    channels = [header['label'] for header in headers]
    scores = pd.read_csv("/data/thmmy/Year_4/Semester_8/PTES/Project/dataset/recordings/SN001_sleepscoring.txt")
    idx = scores[' Annotation'].apply(
        lambda s: re.fullmatch("Sleep.*", s.strip()) is not None
    )
    scores = scores[idx]
    epochs = [
        signals[:, i*dt*fs : (i+1)*dt*fs] for i,_ in enumerate(scores[' Recording onset'])
    ]
    stages = [
        s.split(" ")[-1] for s in scores[' Annotation']
    ]
    return epochs, stages


