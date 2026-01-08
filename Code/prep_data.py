import os
from Bio import SeqIO
import pandas as pd
import glob
import argparse
# set root of the python file
ROOT = os.path.dirname(os.path.abspath(__file__))

# path for folders of interest
folder_M1 = "E1402-Plaque2-M1-20250325"
folder_M2 = "E1402-Plaque2-M2-20250325"


# Raw data
def raw_data(PATH, folder_path, folder_name, overwrite=True):


    result_path = os.path.join(PATH, f"{folder_name}_raw.csv")
    # Optionnel: repartir de zéro pour éviter une ancienne colonne d'index
    if overwrite and os.path.exists(result_path):
        os.remove(result_path)
    first_write = not os.path.exists(result_path)

    for file_path in glob.glob(os.path.join(folder_path, "*.fsa")):
        record = SeqIO.read(file_path, "abi")
        filename = os.path.basename(file_path)
        print("Read :", filename)

        raw = record.annotations["abif_raw"]

        # Check for data that is 4961 long
        signal_keys = [k for k, v in raw.items() if isinstance(v, tuple) and len(v) == 4961]
        # print(signal_keys)


        channels = {

            "channel_1": record.annotations["abif_raw"][signal_keys[0]],
            "channel_2": record.annotations["abif_raw"][signal_keys[1]],
            "channel_3": record.annotations["abif_raw"][signal_keys[2]],
            "channel_4": record.annotations["abif_raw"][signal_keys[3]],
            "channels_5": record.annotations["abif_raw"][signal_keys[4]]
        }

        raw_df = pd.DataFrame(channels)
        raw_df.insert(0, "filename", filename)
        raw_df.to_csv(
            result_path,
            sep=";",
            header=first_write,
            mode="a",
            index=False,
        )
        if first_write:
            first_write = False
# raw_data(ROOT, folder_M2, "M2", overwrite=True)


# Raw data (df) into a Raw data (dictionnary of df)
def prep_tsfresh(
    ROOT,
    filename,
    write_tsfresh=True,
    tsfresh_filename=None,
    save_raw_dict=False,
):
    """
    Lit le CSV en chunks et retourne un dictionnaire `raw_dict` mappant
    `filename` -> DataFrame (colonnes: channels).

    Si `write_tsfresh` est True, convertit aussi en format long compatible tsfresh
    en flux (streaming) et écrit progressivement dans un fichier CSV sans charger
    l'ensemble en mémoire.

    Si `save_raw_dict` est True, sérialise également le dictionnaire `raw_dict`
    en un fichier pickle (par défaut: `<base>_raw_dict.pkl` à côté des CSV).

    Retourne: (raw_dict, tsfresh_csv_path | None)
    """
    DATA = os.path.join(ROOT, "..", "Data")
    DATA = os.path.normpath(DATA)

    result_path = os.path.join(DATA, filename)

    # Fichier de sortie tsfresh (séparé pour ne pas écraser la source)
    tsfresh_path = None
    if write_tsfresh:
        if tsfresh_filename is None:
            base, ext = os.path.splitext(filename)
            tsfresh_filename = f"{base}_tsfresh.csv"
        tsfresh_path = os.path.join(DATA, tsfresh_filename)
        # Réinitialiser le fichier s'il existe déjà
        if os.path.exists(tsfresh_path):
            os.remove(tsfresh_path)

    # Préparer le chemin de sauvegarde pour le dictionnaire brut si demandé
    raw_dict_path = None
    if save_raw_dict:
        base, ext = os.path.splitext(filename)
        raw_dict_filename = f"{base}_dict.pkl"
        raw_dict_path = os.path.join(DATA, raw_dict_filename)
        if os.path.exists(raw_dict_path):
            os.remove(raw_dict_path)

    raw_dict = {}
    # Suivi de l'offset temporel par fichier pour un index "time" continu
    time_offsets = {}
    header_written = False

    for chunk in pd.read_csv(result_path, sep=";", chunksize=5000):
        for name, df_sub in chunk.groupby("filename"):
            df = df_sub.drop(columns=["filename"]).reset_index(drop=True)
            # Accumulation en mémoire du brut par fichier (si souhaité ensuite)
            if name not in raw_dict:
                raw_dict[name] = df.copy()
            else:
                raw_dict[name] = pd.concat([raw_dict[name], df], ignore_index=True)
            
            # Ecriture en flux du format long tsfresh
            if write_tsfresh:
                start = time_offsets.get(name, 0)
                n = len(df)
                if n > 0:
                    time_col = pd.Series(range(start, start + n), name="time")
                    time_offsets[name] = start + n
                    long_df = (
                        df.assign(time=time_col)
                          .melt(id_vars=["time"], var_name="channel", value_name="value")
                    )
                    long_df.insert(0, "id", name + "_" + long_df["channel"].astype(str))
                    long_df = long_df[["id", "time", "value"]]
                    long_df.to_csv(
                        tsfresh_path,
                        sep=";",
                        index=False,
                        header=not header_written,
                        mode="a",
                    )
                    header_written = True

    # Sauvegarde du dictionnaire brut sur disque si demandé
    if save_raw_dict:
        pd.to_pickle(raw_dict, raw_dict_path)

    return raw_dict, tsfresh_path

#filename = "M1_raw.csv"
#prep_tsfresh(ROOT,filename,write_tsfresh=False,tsfresh_filename=None,save_raw_dict=True)

def build_labels_stream(
        ROOT,
        files,
        output_filename,
        raw_folder,
        chunksize=5000,
        overwrite=True
):
    """
    Écrit un fichier de sortie agrégé (colonnes: plant_number, channel, Y1, Y2)
    en mode streaming par chunks pour éviter la saturation mémoire.

    Paramètres:
    - ROOT: racine du script
    - files: liste de chemins relatifs sous Data
    - output_filename: nom du fichier CSV de sortie sous Data
    - chunksize: taille des chunks pour la lecture pandas
    - overwrite: si True, supprime le fichier de sortie existant
        - raw_folder: dossier des fichiers bruts (.fsa) pour mapper le rang (plant_number)
            vers le nom de fichier correspondant, selon l'ordre des fichiers du dossier.

        Retourne: chemin absolu du fichier de sortie.

        Remarque:
        - Remplace la colonne "plant" par "filename" qui contient le nom
            du fichier .fsa correspondant (ordre déterminé par le tri alphanumérique
            des fichiers du dossier `raw_folder`).
    """
    DATA = os.path.join(ROOT, "..", "Data/label")
    DATA = os.path.normpath(DATA)
    output_path = os.path.join(DATA, output_filename)

    if overwrite and os.path.exists(output_path):
        os.remove(output_path)

    # Construire le mapping rang -> nom de fichier .fsa depuis le dossier brut
    raw_dir = os.path.normpath(os.path.join(ROOT, "..", "Data", raw_folder))
    fsa_paths = sorted(glob.glob(os.path.join(raw_dir, "*.fsa")))
    fsa_names = [os.path.basename(p) for p in fsa_paths]
    # Mapping 1-based (1 -> premier fichier, 2 -> deuxième, ...)
    index_to_filename = {i + 1: name for i, name in enumerate(fsa_names)}

    header_written = False
    for rel_path in files:
        file_path = os.path.join(DATA, rel_path)
        channel_name = os.path.splitext(os.path.basename(file_path))[0]

        for chunk in pd.read_csv(
            file_path,
            sep=";",
            engine="python",
            quotechar='"',
            na_values=["NA"],
            on_bad_lines="skip",
            chunksize=chunksize,
        ):
            # Renommer et ajouter la colonne channel
            chunk = chunk.rename(columns={
                "markA.1": "Y1",
                "markA.2": "Y2",
            })
            chunk["channel"] = channel_name

            # Construire directement la sortie avec filename mappé (en une fois)
            out = pd.DataFrame({
                "filename": pd.to_numeric(chunk["plant"], errors="coerce").astype("Int64").map(index_to_filename),
                "channel": chunk["channel"],
                "Y1": pd.to_numeric(chunk["Y1"], errors="coerce"),
                "Y2": pd.to_numeric(chunk["Y2"], errors="coerce"),
            })

            out.to_csv(
                output_path,
                sep=";",
                index=False,
                header=not header_written,
                mode="a",
            )
            header_written = True

    return output_path

#files = ["M2/channel_1.csv", "M2/channel_2.csv", "M2/channel_3.csv", "M2/channel_4.csv"]
#build_labels_stream(ROOT, files, output_filename="M2_labels_stream.csv", raw_folder= folder_M2)

"""
def split_dataset(ROOT,tsfresh_df, label, merge=True):
    if merge:
        
    X_train, X_test, y_train, y_test = train_test_split(data, label, random_state=20240707)
def get_signal(raw_dict, filename, channel):
    if filename not in raw_dict:
        raise KeyError(f"Unknown filename: {filename}")
    df = raw_dict[filename]
    if channel not in df.columns:
        raise KeyError(f"Unknown channel: {channel}")
    return df[channel].reset_index(drop=True)

# Feature Extraction
# def feature_extraction(raw_dict):

# associer X et Y

# Main
# def main()
    DATA = os.path.join(ROOT, "..", "Data")
    DATA = os.path.normpath(DATA)
    PATH = os.path.join(DATA, folder_path)
# argpase
"""