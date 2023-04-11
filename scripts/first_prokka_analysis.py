import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

PATH = "//wsl.localhost/Ubuntu/home/vitt/MICROBIAL_GENOMICS/project/data/files"

hypothetical_dict = {}
CDS_dict = {}
hp_CDS_dict = {}

for file in os.listdir(PATH):
    MAG_dir_path = os.path.join(PATH, file)
    short_path = file

    # checking if it is a dir
    if os.path.isdir(MAG_dir_path):

        for file in os.listdir(MAG_dir_path):
            tsv_file = os.path.join(MAG_dir_path, file) 

            # find all TSV files
            if tsv_file.endswith(".tsv"):

                print(tsv_file)
                
                tsv_file = pd.read_table(tsv_file, index_col=0)
                CDS_dict[short_path] = tsv_file.ftype.value_counts().CDS
                temp_dict = dict(tsv_file["product"].value_counts())
                hypothetical_dict[short_path] = temp_dict["hypothetical protein"]

                counter = 0
                for index, row in tsv_file.iterrows():
                  if row["product"] != "hypothetical protein" and row["ftype"] == "CDS":
                    counter += 1
                hp_CDS_dict[short_path] = counter

# print("\n", CDS_dict)
# print("\n\n\n", hypothetical_dict)

sum_CDS = 0
for key in CDS_dict:
    sum_CDS += CDS_dict[key]

sum_hp = 0
for key in hypothetical_dict:
    sum_hp += hypothetical_dict[key]

sum_hp_CDS = 0
for key in hp_CDS_dict:
    sum_hp_CDS += hp_CDS_dict[key]

# print("avg of CDS:", sum_CDS/len(CDS_dict))
# print("avg of hypothetical protein:", sum_hp/len(hypothetical_dict))

plt.plot(list(CDS_dict.keys()), list(CDS_dict.values()), "-", 
         list(CDS_dict.keys()), list(hypothetical_dict.values()), "-",
         list(CDS_dict.keys()), list(hp_CDS_dict.values()), "-")
plt.xticks([])
plt.gca().legend(('CDS','HP', 'KP'))
plt.title("prokka summary for SGB4964")
plt.axhline(y=sum_CDS/len(CDS_dict), color='tab:blue', alpha=0.5)
plt.axhline(y=sum_hp/len(hypothetical_dict), color='tab:orange', alpha=0.5)
plt.axhline(y=sum_hp_CDS/len(hp_CDS_dict), color='tab:green', alpha=0.5)
plt.text(0, sum_CDS/len(CDS_dict) + 150,'avg =' + str(round(sum_CDS/len(CDS_dict))), color='tab:blue')
plt.text(1, sum_hp/len(hypothetical_dict) - 170,'avg =' + str(round(sum_hp/len(hypothetical_dict))), color='tab:orange')
plt.text(0, sum_hp_CDS/len(hp_CDS_dict) + 150,'avg =' + str(round(sum_hp_CDS/len(hp_CDS_dict))), color='tab:green')
plt.show()

