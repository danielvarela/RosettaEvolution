
# script for create the paper violin plots and concatenate them in different groups as it should go in the paper

import os
import shutil

def make_folder(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def make_group(name, group):
    make_folder(name)
    positions = "a b c d".split()
    counter = 0
    for p in group:
        energy_name = "violin_plot_Supp_energy_"+p+".png"
        shutil.copy(energy_name, "./"+name+"/"+positions[counter]+"_"+energy_name)
        counter = counter + 1
        rmsd_name = "violin_plot_Supp_rmsd_"+p+".png"
        shutil.copy(rmsd_name, "./"+name+"/"+positions[counter]+"_"+rmsd_name)
        counter = counter + 1
                   
group_1 = ["1kpeA", "1rnbA"]
group_2 = ["1dhn", "1eyvA"]
group_3 = ["1cg5B", "1lis"]

os.system("python3 pyscripts/finder.py paper_bilevel_database_aux.csv --graph cluster --out_folder rep_again_3colores --prot_set all --violin_type rmsd")
os.system("python3 pyscripts/finder.py paper_bilevel_database_aux.csv --graph cluster --out_folder rep_again_3colores --prot_set all --violin_type energy ")



make_group("group_1", group_1)
make_group("group_2", group_2)
make_group("group_3", group_3)

os.system("montage -tile 2x2 group_1/*.png -mode concatenate -geometry -8 NEW_paper_violin_combined1.png")

os.system("montage -tile 2x2 group_2/*.png -mode concatenate -geometry -8 NEW_paper_violin_combined2.png")
os.system("montage -tile 2x2 group_3/*.png -mode concatenate -geometry -8 NEW_paper_violin_combined3.png")
