import os
from Bio import SeqIO
import numpy as np
import pandas as pd
import glob
import argparse
#set root of the python file
ROOT = os.path.dirname(os.path.abspath(__file__))

#path for folders of interest
folder_M1 = "E1402-Plaque2-M1-20250325"
folder_M2 = "E1402-Plaque2-M2-20250325"

##Raw data
def raw_data(ROOT, folder_path, folder_name):

    DATA = os.path.join(ROOT, "..", "Data")
    DATA = os.path.normpath(DATA)
    path = os.path.join(DATA,folder_path)

    result_path = os.path.join(DATA,f"{folder_name}_raw.csv")
    first_write = not os.path.exists(result_path)

    for file_path in glob.glob(os.path.join(path, "*.fsa")):
        record = SeqIO.read(file_path, "abi")
        filename=os.path.basename(file_path)
        print("Read :", filename)

        raw = record.annotations["abif_raw"]

        #Check for data that is 4961 long
        signal_keys = [k for k, v in raw.items() if isinstance(v, tuple) and len(v) == 4961]
        #print(signal_keys)


        channels = {
            "channel_1": record.annotations["abif_raw"][signal_keys[0]],
            "channel_2": record.annotations["abif_raw"][signal_keys[1]],
            "channel_3": record.annotations["abif_raw"][signal_keys[2]],
            "channel_4": record.annotations["abif_raw"][signal_keys[3]],
            "channels_5": record.annotations["abif_raw"][signal_keys[4]]
        }

        raw_df = pd.DataFrame(channels)
        raw_df.insert(0, "filename", filename)
        raw_df.to_csv(result_path, sep=";", header=first_write, mode="a")
        if first_write:
            first_write=False
    
raw_data(ROOT, folder_M2, "M2")

## Raw data (df) into a Raw data (dictionnary of df)
def prep_tsfresh(DATA):
    result_path = os.path.join(DATA, "raw_data.csv")

    raw_dict = {}

    for chunk in pd.read_csv(result_path, sep=";", chunksize=5000):
        for name, df_sub in chunk.groupby("filename"):
            if name not in raw_dict:
                raw_dict[name] = df_sub.drop(columns=["filename"])
            else:
                raw_dict[name] = pd.concat([raw_dict[name], df_sub.drop(columns=["filename"])], ignore_index=True)
    return raw_dict
##Feature Extraction

##Main
def main():
##argpase