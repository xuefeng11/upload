import os
import pandas as pd
import numpy as np
import pybel

folder_path_2="."

for file in os.listdir(folder_path_2):
    file_path = folder_path_2 + "/" + file
    if "DS_Store" in file_path:
        continue
    try:
        if file.endswith(".can"):
            df = pd.read_csv(file_path, sep='\t', names=["canonical_smile", "ID"])
            print("processing! ",file_path)

            print(file," info ", df.info())
        else:
            continue

        df = df.drop_duplicates("canonical_smile")
        df.to_csv("unique_"+file,sep="\t")
        print("unique ",file,df.info())

        print("#####################################")

    except Exception as e:
        print("Error!", e)
        exit()

