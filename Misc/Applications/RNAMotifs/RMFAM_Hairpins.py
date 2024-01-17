from Bio import AlignIO
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

def load_internals(path_sr1,path_sr2,path_kturn1,path_kturn2,path_cloop): #I started working on a version for Internal Loops aswell, after thorough investigation of their stockholm files I came to the realization that this is fruitless as their structures are to inhomogenous to efficiently program this (its faster for me to just read through the files and write down the sequences.)
    sr1_file=load_stk(path_sr1)
    sr2_file=load_stk(path_sr2)
    kturn1_file=load_stk(path_kturn1)
    kturn2_file=load_stk(path_kturn2)
    cloop_file=load_stk(path_cloop)
    return sr1_file,sr2_file,kturn1_file,kturn2_file,cloop_file

def write_sequences(loops_set,single_letter_abbreviation,path_rmfam_hairpins):
    with open(path_rmfam_hairpins,"a") as csv_file:
        for entry in loops_set:
            entry_abb=entry+"+{abb}\n".format(abb=single_letter_abbreviation)
            csv_file.write(entry_abb)

if __name__ == "__main__":
    print("REMEMBER TO DELETE THE OLD HairpinRMFAM.csv FILE WHENEVER YOU RUN THIS PROGRAM AS THIS JUST APPENDS TO THE EXISTING ONE OTHERSWISE! WHILE THIS DONT MAKE A DIFFERENCE FOR PREDICTIONS IF THERE ARE CHANGES IN THE STOCKHOLM FILES YOU WILL KEEP THE OLD SEQUENCES IN THERE!")
    tloops_stk,gnra_stk,uncg_stk=load_hairpins("Misc/Applications/RNAMotifs/Data/RM00024.stockholm.txt","Misc/Applications/RNAMotifs/Data/RM00008.stockholm.txt","Misc/Applications/RNAMotifs/Data/RM00029.stockholm.txt")
    tloops=get_sequences_hairpins(tloops_stk,21,30)
    gnras=get_sequences_hairpins(gnra_stk,35,39)
    uncgs=get_sequences_hairpins(uncg_stk,24,28)
    write_sequences(tloops,"T","Misc/Applications/RNAMotifs/Loops/HairpinRMFAM.csv")
    write_sequences(gnras,"G","Misc/Applications/RNAMotifs/Loops/HairpinRMFAM.csv")
    write_sequences(uncgs,"U","Misc/Applications/RNAMotifs/Loops/HairpinRMFAM.csv")