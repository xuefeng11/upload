import os
import pandas as pd
import numpy as np
import pybel


path="DUDE"
smile_arr=list()
id_arr = list()
for folder in os.listdir(path):
    folder_path=path+"/"+folder
    if "DS_Store" in folder:
        continue
    for dir in os.listdir(folder_path):
        folder_path_2=folder_path+"/"+dir
        if "DS_Store" in folder_path_2:
            continue
        for file in os.listdir(folder_path_2):
            file_path = folder_path_2 + "/" + file
            if "DS_Store" in file_path:
                continue

            try:
                if file.endswith("actives_final.ism"):
                    df = pd.read_csv(file_path, sep=' ',names=["smile","num","ID"])
                elif file.endswith("decoys_final.ism"):
                    df = pd.read_csv(file_path, sep=' ', names=["smile","ID"])
                else:
                    continue

                print(file_path)
                smile_arr = smile_arr + df["smile"].tolist()
                id_arr = id_arr + df["ID"].tolist()

            except Exception as e:
                print("Error!",e)
                exit()


data = {"smile":smile_arr,"name": id_arr}
df_out = pd.DataFrame(data)
print("df raw:",df_out.info())
df_out.to_csv("DUDE_RAW_SMILE.csv",index=False)
print("---------------------------")


df_out=df_out.drop_duplicates("smile")
print("df unique smile:",df_out.info())
df_out.to_csv("DUDE_UNIQUE_SMILE.csv", index=False)
print("---------------------------")

drug=df_out["smile"]
canonical_smiles = [pybel.readstring("smi", smile).write("can").rstrip() for smile in drug]
data = {"canonical_smile":canonical_smiles,"ID": df_out['name']}
df_out = pd.DataFrame(data)
print("df canonical smile:",df_out.info())
print("---------------------------")

df_out=df_out.drop_duplicates("canonical_smile")
print("df unique canonical smile:",df_out.info())
df_out.to_csv("DUDE_UNIQUE_canonical.csv", index=False)
print("---------------------------")
