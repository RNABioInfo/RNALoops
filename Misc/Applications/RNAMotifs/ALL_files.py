import os.path
import inspect
import glob
import re
from collections import Counter

#THIS ALGORITHM TAKES THE PREWRITTEN/GENERATED .csv FILES IN BULGES HAIRPINS AND INTERNALS AND WRITES THE ALL FILES FROM THEM (COMBINES THEM) FIRST DELETES EXISTING ALL FILES, THEN REWRITES
#THEM FROM THE EXISTING BGSU AND RMFAM FILES. IF YOU WANT TO ADD MORE DATABASES THERE ALWAYS NEEDS TO BE A FORWARD REVERSE AND BOTH FILES. THIS CAN JUST ADD THESE INTO THE ALL FILES
#ALSO:
#with check for dupes this algorithms also contains a check to see whether any sequences have two associated motifs, as of today (23.01.2024, BGSU Version 3.77) there is
# only one sequence: GUGA (G/U), I'll take the U version out in some way to avoid unpredictable behavior in the motif prediction --> Is taken out of BGSU_forward, BGSU_reverse and BGSU_both.

def get_filepaths():
    filename  = inspect.getframeinfo(inspect.currentframe()).filename
    path      = os.path.dirname(os.path.abspath(filename))
    b_loops_path_B= path + "/Loops/Bulges/BGSU*"
    i_loops_path_B= path + "/Loops/Internals/BGSU*"
    h_loops_path_B= path + "/Loops/Hairpins/BGSU*"
    b_loops_path_R= path + "/Loops/Bulges/RMFAM*"
    i_loops_path_R= path + "/Loops/Internals/RMFAM*"
    h_loops_path_R= path + "/Loops/Hairpins/RMFAM*"
    b_files=glob.glob(b_loops_path_B) + glob.glob(b_loops_path_R)
    h_files=glob.glob(h_loops_path_B) + glob.glob(h_loops_path_R)
    i_files=glob.glob(i_loops_path_B) + glob.glob(i_loops_path_R)
    return b_files,h_files,i_files,path
    
def find_sets(file_path_list):
    forward = re.compile("_forward.csv")
    reverse = re.compile("_reverse.csv")
    both    = re.compile("_both.csv")
    fw_paths=list(filter(forward.search,file_path_list))
    re_paths=list(filter(reverse.search,file_path_list))
    bo_paths=list(filter(both.search,file_path_list))
    return fw_paths,re_paths,bo_paths

def combine_to_ALL(source):
    motif_list=[]
    for entry in source:
        with open(entry,"r") as source_file:
            list_source_file = source_file.readlines()
            motif_list.extend(list_source_file)
    return motif_list

def write_to_csv(source_list,cwd,path):
    dupe_check(source_list)
    file_path2=cwd+path
    if os.path.isfile(file_path2) == True:
        os.remove(file_path2)
    else:pass
    with open(file_path2,"w") as output_file:
        for entry in source_list:
            output_file.write(entry)

def dupe_check(check_list):
    dic={}
    for entry in check_list:
        mots=[]
        split=entry.split("+")
        mots.append(split[1])
        if split[0] not in dic.keys():
            dic[split[0]]=mots
        else:
            dic[split[0]].extend(mots)
    for entry in dic:
        if len(set(dic[entry])) > 1:
            print(entry,set(dic[entry]))

if __name__ == "__main__":
    bulges,hairpins,internals,path=get_filepaths()
    forward_b,reverse_b,both_b=find_sets(bulges)
    forward_h,reverse_h,both_h=find_sets(hairpins)
    forward_i,reverse_i,both_i=find_sets(internals)
    ALL_fw_b=combine_to_ALL(forward_b)
    write_to_csv(ALL_fw_b,path,"/Loops/Bulges/ALL_forward.csv")
    ALL_re_b=combine_to_ALL(reverse_b)
    write_to_csv(ALL_re_b,path,"/Loops/Bulges/ALL_reverse.csv")
    ALL_bo_b=combine_to_ALL(both_b)
    write_to_csv(ALL_bo_b,path,"/Loops/Bulges/ALL_both.csv")
    ALL_fw_h=combine_to_ALL(forward_h)
    write_to_csv(ALL_fw_h,path,"/Loops/Hairpins/ALL_forward.csv")
    ALL_re_h=combine_to_ALL(reverse_h)
    write_to_csv(ALL_re_h,path,"/Loops/Hairpins/ALL_reverse.csv")
    ALL_bo_h=combine_to_ALL(both_h)
    write_to_csv(ALL_bo_h,path,"/Loops/Hairpins/ALL_both.csv")
    ALL_fw_i=combine_to_ALL(forward_i)
    write_to_csv(ALL_fw_i,path,"/Loops/Internals/ALL_forward.csv")
    ALL_re_i=combine_to_ALL(reverse_i)
    write_to_csv(ALL_re_i,path,"/Loops/Internals/ALL_reverse.csv")
    ALL_bo_i=combine_to_ALL(both_i)
    write_to_csv(ALL_bo_i,path,"/Loops/Internals/ALL_both.csv")