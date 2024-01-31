from Bio import AlignIO
import os.path

def load_stk(path_loops):
    file=AlignIO.parse(path_loops,format="stockholm")
    return file
         
def get_sequences_hairpins(loop_file,closing_base1,closing_base2):
    seq_list=[]
    for entry in loop_file:
        for entry2 in entry._records:
            seq=list(entry2.seq[closing_base1:closing_base2])
            List_seq_stripped=[x.upper() for x in seq if x != "-"]
            if "N" not in List_seq_stripped:
                string_seq=''.join(List_seq_stripped)
                seq_list.append(string_seq)
            else:
                pass
    return set(seq_list)

def load_hairpins(path_tloops,path_gnras,path_uncgs):
    tloop_file=load_stk(path_tloops)
    gnra_file=load_stk(path_gnras)
    uncg_file=load_stk(path_uncgs)
    return tloop_file,gnra_file,uncg_file

def write_sequences(loops_set,single_letter_abbreviation,path_rmfam_hairpins_forward,path_rmfam_hairpins_reverse,path_rmfam_hairpins_both):
    for entry in loops_set:
        entry_abb=entry+"+{abb}\n".format(abb=single_letter_abbreviation)
        entry_reverse_abb=entry[::-1]+"+{abb}\n".format(abb=single_letter_abbreviation)
        with open(path_rmfam_hairpins_forward,"a") as csv_file:
            csv_file.write(entry_abb) 
        with open(path_rmfam_hairpins_reverse,"a") as csv_file:
            csv_file.write(entry_reverse_abb)    
        with open(path_rmfam_hairpins_both,"a") as csv_file:
            csv_file.write(entry_abb)
            csv_file.write(entry_reverse_abb)
    return 0

def load_internals(path_sr1,path_sr2,path_kturn1,path_kturn2,path_cloop): #I started working on a version for Internal Loops aswell, after thorough investigation of their stockholm files I came to the realization that this is fruitless as their structures are to inhomogenous to efficiently program this (its faster for me to just read through the files and write down the sequences.)
    sr1_file=load_stk(path_sr1)
    sr2_file=load_stk(path_sr2)
    kturn1_file=load_stk(path_kturn1)
    kturn2_file=load_stk(path_kturn2)
    cloop_file=load_stk(path_cloop)
    return sr1_file,sr2_file,kturn1_file,kturn2_file,cloop_file

def get_sequences_sr1(sr1_file):#Yeah i kept trying it with the Internals but this really aint working, might need to do it for paper release anyways but the heterogenous nature of the .stk files makes it really impractical.
    for entry in sr1_file:
        motpos=entry._per_col_annotations['secondary_structure'].index("[")
        motpos2=entry._per_col_annotations['secondary_structure'].index("]")
        for seq in entry._records:
            print(seq.seq[motpos:motpos+4]+"$"+seq.seq[motpos2:motpos2+5])
        print(entry._per_col_annotations['secondary_structure'])

def reverse_internals(path_to_forward_internals,path_to_reverse_internals, path_to_both_internals):
    if os.path.isfile(path_to_reverse_internals) == True:
        os.remove(path_to_reverse_internals)
        os.remove(path_to_both_internals)
    with open(path_to_forward_internals,"r") as input_file:
        for line in input_file:
            split_line=line.split("+")
            seq=split_line[0]
            reverse_seq=seq[::-1]
            reverse_seq_abb=reverse_seq+"+{abb}".format(abb=split_line[1])
            with open(path_to_reverse_internals,"a") as output_file1:
                    output_file1.write(reverse_seq_abb)
            with open(path_to_both_internals,"a") as output_file2:
                    output_file2.write(line)
                    output_file2.write(reverse_seq_abb) 
    return 0

def remove_old_files(hl_f,hl_r,hl_b):
    if os.path.isfile(hl_f) == True:
        os.remove(hl_f)
        os.remove(hl_r)
        os.remove(hl_b)

if __name__ == "__main__":
    #sr1,sr2,kturn1,kturn2,cloop=load_internals("Misc/Applications/RNAMotifs/Data/RM00018.stockholm.txt","Misc/Applications/RNAMotifs/Data/RM00019.stockholm.txt","Misc/Applications/RNAMotifs/Data/RM00010.stockholm.txt","Misc/Applications/RNAMotifs/Data/RM00011.stockholm.txt","Misc/Applications/RNAMotifs/Data/RM00003.stockholm.txt")
    reverse_internals("/home/ubuntu/RNALoops/Misc/Applications/RNAMotifs/Loops/Internals/RMFAM_forward.csv","/home/ubuntu/RNALoops/Misc/Applications/RNAMotifs/Loops/Internals/RMFAM_reverse.csv","/home/ubuntu/RNALoops/Misc/Applications/RNAMotifs/Loops/Internals/RMFAM_both.csv")
    tloops_stk,gnra_stk,uncg_stk=load_hairpins("Misc/Applications/RNAMotifs/Data/RM00024.stockholm.txt","Misc/Applications/RNAMotifs/Data/RM00008.stockholm.txt","Misc/Applications/RNAMotifs/Data/RM00029.stockholm.txt")
    tloops=get_sequences_hairpins(tloops_stk,21,30)
    gnras=get_sequences_hairpins(gnra_stk,35,39)
    uncgs=get_sequences_hairpins(uncg_stk,24,28)
    remove_old_files("Misc/Applications/RNAMotifs/Loops/Hairpins/HairpinRMFAM_forward.csv", "Misc/Applications/RNAMotifs/Loops/Hairpins/HairpinRMFAM_reverse.csv", "Misc/Applications/RNAMotifs/Loops/Hairpins/HairpinRMFAM_both.csv")
    write_sequences(tloops,"T" ,"Misc/Applications/RNAMotifs/Loops/Hairpins/HairpinRMFAM_forward.csv", "Misc/Applications/RNAMotifs/Loops/Hairpins/HairpinRMFAM_reverse.csv", "Misc/Applications/RNAMotifs/Loops/Hairpins/HairpinRMFAM_both.csv")
    write_sequences(gnras ,"G" ,"Misc/Applications/RNAMotifs/Loops/Hairpins/HairpinRMFAM_forward.csv", "Misc/Applications/RNAMotifs/Loops/Hairpins/HairpinRMFAM_reverse.csv", "Misc/Applications/RNAMotifs/Loops/Hairpins/HairpinRMFAM_both.csv")
    write_sequences(uncgs ,"U" ,"Misc/Applications/RNAMotifs/Loops/Hairpins/HairpinRMFAM_forward.csv", "Misc/Applications/RNAMotifs/Loops/Hairpins/HairpinRMFAM_reverse.csv", "Misc/Applications/RNAMotifs/Loops/Hairpins/HairpinRMFAM_both.csv")