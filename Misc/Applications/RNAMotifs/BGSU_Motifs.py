###################### READ BEFORE STARTING THE PROGRAM##############################
#
#This program was written to work in the RNALoops folder downloadable from my GitHub page
#MSebeke. File paths are made in a way that assumes that the program is located in the
#RNALoops folder but can be ran from anywhere. The important part is that it is located
#in a folder with two other folders: A ./Data/ folder that contains three different files:
#1. hl_x.xx.json 2. il_x.xx.json 3. Loop_List.txt
#1 and 2 are the downloadable .json files from the BGSU database, specifically this program was written
#when the third iteration of the BGSU was around, so later versions might break this based on different structures
# in the .json or a changes in the API returns. The third file, Loop_List.txt, specifies the selected loops that 
#the user wants to have the sequences written in .csv file for usage with our RNALoops.gap.
#The sequences are written to the second important folder, which would be ./Loops, this folder is the place the .csv files
#are written to and for this reason has to exist as I have to still implement a check for whether the folder exists.
#For additional information on the Loop_List.txt, refer to the relevant_loops(path_to_loops) function or reach out to me.
#WRITTEN ON 28.09.2023 BY MARIUS SEBEKE AS PART OF THE RNALOOPS PROJECT, FOR QUESTIONS REACH OUT TO:  marius.sebeke@ibmg.uni-stuttgart.de 
###################################################################################

###### Imports ######
import inspect
import os.path
import glob
import re
from time import sleep
import json
from collections import Counter
import requests
import sys

##### Functions #####

def get_filepaths(): #Finds the paths to the directory where this python file is stored ("something/something/RNALoops/Extensions/RNAMotifs") and then builds the paths to the necessary downloadable json files from the BGSU 
    filename  = inspect.getframeinfo(inspect.currentframe()).filename
    path      = os.path.dirname(os.path.abspath(filename))
    data_path = path+"/Data/*"
    file_list=glob.glob(data_path)
    hairpins_re=re.compile(".*hl_")
    internals_re=re.compile(".*il_")
    loop_list_re=re.compile(".*Loop_List.txt")
    hairpins_path=list(filter(hairpins_re.match,file_list))[0]
    internals_path=list(filter(internals_re.match,file_list))[0]
    loop_list_path=list(filter(loop_list_re.match,file_list))[0]
    return hairpins_path,internals_path,loop_list_path
    
def relevant_loops(path_to_loops): #Takes the path to the Loop List.txt, which is the file that should be saved in "".../RNALoops/Extensions/RNAMotifs/Data/". This file specifies the Motif IDs from the BGSU that are to be used
    relevant_loops_dict={} #in the creation of the motif sequences. This can be modified to only contain certain motifs if only those are of interest. File formatting is easy: MOTIF_ID (\t or 4 spaces) ABBREVIATION (e.g. HL_79299    3)
    with open(path_to_loops,"r") as file:
        Lines=file.readlines()
    for Line in Lines:
        Stripped_Line=Line.strip()
        Split_Line=re.split("    ",Stripped_Line)
        relevant_loops_dict[Split_Line[0]]=Split_Line[1]
    return relevant_loops_dict

def get_loop_ids(path_to_hairpins,path_to_internals): #Reads through the hl.json and il.json files and creates a dictionary containing all Motifs with their ID as key and all the corresponding motif instances in a list as the value
    with open(path_to_hairpins, "r") as file:
        hairpins_json=json.load(file)
    with open(path_to_internals,"r") as file:
        internals_json=json.load(file)
    loops_json=hairpins_json+internals_json
    id_motif_dict={}
    for i in range(len(loops_json)):
        Split_ID=re.split("[.]",loops_json[i]["motif_id"])
        id=Split_ID[0]
        id_motif_dict[id]=list(loops_json[i]["alignment"].keys())
    return id_motif_dict

def API_prep(idm_d,cl_d): #Preps the API calls by creating a dictionary with a Motif ID as key with all it's corresponding motif instance API calls in a list as value
    id_api_dict={}
    for id in idm_d: #All IDs with all their motif instances
        api_list=[]
        if id in cl_d.keys(): #Filter that only ID's that are in the Loop_List pass!
            for motif_instance in idm_d[id]:
                api_request="http://rna.bgsu.edu/correspondence/pairwise_interactions_single?selection_type=loop_id&selection="+motif_instance
                api_list.append(api_request)
        if len(api_list) != 0: #Just in case filter that prevents an empty api_list from being handed off to the rest of the program which would just lead to an error.
            id_api_dict[id]=api_list
    return id_api_dict

def API(ida_d): #Function that actually executes the API calls. Takes the API call list as a request list and removes the entries that successfully executed. If a call fails it gives the server a 2 second break before retrying. 
    id_response_dict={} #Successful calls result in the answer being saved as a list in a dictionary with the ID as key and a list of the Answer lists. 
    for id in ida_d:
        response_list=[]
        request_list=ida_d[id]
        while len(request_list) != 0:
            for api_call in request_list:
                response=requests.get(api_call)
                if response.status_code == 200:
                    request_list.remove(api_call)
                    decoded_response=response.content.decode()
                    response_list.append(decoded_response)
                else:
                    sleep(5)
        id_response_dict[id]=response_list
    return id_response_dict

def API_postprocess_part1(idr): #Postprocessing of API responses. Iterates through the returned lists with two effects: Finds all nucleotides  (except the very first one cause it is preceded by a <\br> which isn't caught by the regex)
    RegEx=re.compile("[a-zA-Z0-9_]+[|]") # And it pops away the very last nucleotide. These removals are both planned as these nucleotides ALWAYS belong to the closing base pair which we do not need for the motif recognition.
    id_nucleotides_dict={} #After this function the Hairpins are already basically done as the closing basepair is removed here. For Internals we need a second step as this removes only one of their two closing base pairs.
    for id in idr: #The second closing basepair is removed in part2.
        list_of_nucleotide_lists=[]
        for returned in idr[id]:
            listed_return = returned.split()
            Nucleotides=list(filter(RegEx.match,listed_return))
            Nucleotides.pop()
            list_of_nucleotide_lists.append(Nucleotides)
        id_nucleotides_dict[id]=list_of_nucleotide_lists
    return id_nucleotides_dict

def API_postprocess_part2(idn): #This function takes in the already worked on id_nucleotide dictionary that was created in part1. It removes the second closing base pair for internal loops as well as filter out whether a sequence belongs
    Hairpin_sequences=[] #To Internal Loop or a Bulge Loop, as both of these are classified as Internals by the BGSU but are treated separately by the RNA secondary structure prediction algorithm.
    Internal_sequences=[] #Separation of the parts for a sequences is done based on the absolute base positions. Sequences are saved in tuples with their single letter abbreviation in lists that sort them into their loop type.
    Bulge_sequences=[]
    for id in idn:
        if re.match("I",id): #Filter out only Internal Loops, this part is for breaking the two parts of the motifs apart and deleting the second closing base pair
            for instance in idn[id]:
                Chainsbool,two_chains_break=internals_chain_check(instance)
                if Chainsbool == False:
                    Break=Find_break(instance) #If there is only one chain the Chainsbool is returned as falls from internals chain check, promting to algorithm to use the Find_Break function for the instance to identify the break
                    FirstPart=instance[0:Break-1]
                    SecondPart=instance[Break+1:]
                else:                          #If the Chainsbool is true the algorithm goes into the else condition and uses the determined chain break from the internals chain check. 
                    FirstPart=instance[0:two_chains_break-1]
                    SecondPart=instance[two_chains_break+1:]
                Sequence,Bulge_bool=Building_sequences_ib_il(FirstPart,SecondPart) #Here we differentiate the Internal Loops from the Bulge Loops. This is done based on the length of the two parts (see function for more detail)
                Result_tuple=(Sequence,id)
                if Bulge_bool == True:
                    Bulge_sequences.append(Result_tuple)
                else:
                    Internal_sequences.append(Result_tuple)
        else:
            for instance in idn[id]:
                Sequence=Building_sequences_base(instance)
                Result_tuple=(Sequence,id)
                Hairpin_sequences.append(Result_tuple)
    return Hairpin_sequences,Internal_sequences,Bulge_sequences

def Writing_sequences(hl,il,bl,cl_d,Rev): #Takes in the tuple lists for all loops, builds lists of strings in the Final formatting. Before writing checks if the files already exist and deletes them if they do. 
    hl_seqs,hl_revseqs=Output_prep(hl,cl_d) #Finally writes all the sequences into the specified files (this isn't split yet since I haven't found a convenient way to not hardcode the paths and this way they're at least all together)
    hl_bothseqs=hl_seqs.union(hl_revseqs)
    il_seqs,il_revseqs=Output_prep(il,cl_d)
    il_bothseqs=il_seqs.union(il_revseqs)
    bl_seqs,bl_revseqs=Output_prep(bl,cl_d)
    bl_bothseqs=bl_seqs.union(bl_revseqs)
    filename  = inspect.getframeinfo(inspect.currentframe()).filename
    path      = os.path.dirname(os.path.abspath(filename))
    if Rev == "1":
        hl_path = path + "/Loops/Hairpins/BGSU_forward.csv"
        il_path = path + "/Loops/Internals/BGSU_forward.csv"
        bl_path = path + "/Loops/Bulges/BGSU_forward.csv"
        remove_previous_files(hl_path,il_path,bl_path)
        write_new_files(hl_path, il_path, bl_path, hl_seqs, il_seqs, bl_seqs)     
    if Rev == "2":
        hl_path = path + "/Loops/Hairpins/BGSU_reverse.csv"
        il_path = path + "/Loops/Internals/BGSU_reverse.csv"
        bl_path = path + "/Loops/Bulges/BGSU_reverse.csv"
        remove_previous_files(hl_path,il_path,bl_path)
        write_new_files(hl_path, il_path, bl_path, hl_revseqs, il_revseqs, bl_revseqs)       
    if Rev == "3":
        hl_path = path + "Loops/Hairpins/BGSU_both.csv"
        il_path = path + "Loops/Internals/BGSU_both.csv"
        bl_path = path + "Loops/Bulges/BGSU_both.csv"
        remove_previous_files(hl_path,il_path,bl_path)
        write_new_files(hl_path, il_path, bl_path, hl_bothseqs, il_bothseqs, bl_bothseqs)      
    if Rev == "4":
        hl_path  = path + "/Loops/Hairpins/BGSU_forward.csv"
        hl_path2 = path + "/Loops/Hairpins/BGSU_reverse.csv"
        hl_path3 = path + "/Loops/Hairpins/BGSU_both.csv"      
        il_path  = path + "/Loops/Internals/BGSU_forward.csv"
        il_path2 = path + "/Loops/Internals/BGSU_reverse.csv"
        il_path3 = path + "/Loops/Internals/BGSU_both.csv"
        bl_path  = path + "/Loops/Bulges/BGSU_forward.csv"
        bl_path2 = path + "/Loops/Bulges/BGSU_reverse.csv"
        bl_path3 = path + "/Loops/Bulges/BGSU_both.csv"
        remove_previous_files(hl_path,il_path,bl_path)
        remove_previous_files(hl_path2,il_path2,bl_path2)
        remove_previous_files(hl_path3,il_path3,bl_path3)
        write_new_files(hl_path ,il_path ,bl_path ,hl_seqs    ,il_seqs    ,bl_seqs)
        write_new_files(hl_path2,il_path2,bl_path2,hl_revseqs ,il_revseqs ,bl_revseqs)
        write_new_files(hl_path3,il_path3,bl_path3,hl_bothseqs,il_bothseqs,bl_bothseqs)
    else: print("Please choose a Reverse Parameter of 1 through 4, see Rev param description for more information")
    return 0

####### Sub-functions #######

def remove_previous_files(hl_p,il_p,bl_p):
    if os.path.isfile(hl_p) == True:
        os.remove(hl_p)
        os.remove(il_p)
        os.remove(bl_p)
    else:pass
    
def write_new_files(hl_p,il_p,bl_p,hl_s,il_s,bl_s):
    with open(hl_p,"x") as file:
        for entry in hl_s:
            file.write(entry+"\n")
    with open(il_p,"x") as file:
        for entry in il_s:
            file.write(entry+"\n")
    with open(bl_p,"x") as file:
        for entry in bl_s:
            file.write(entry+"\n")
             
def internals_chain_check(listed_instance): #Takes a list of Nucleotides, checks if the two Parts of the Internal Loop are parts of different Chains and if they are finds the break based on the chains.
    checklist=[]
    for nucleotide in listed_instance:
        split_n=nucleotide.split("|")
        chainid=split_n[2]
        checklist.append(chainid)
    first_chain_id=checklist[0]
    Chain_counter=Counter(checklist)
    Two_Chains_Break=Chain_counter[first_chain_id] #If there are two Chains this sets the break based on how often the Chain that the first Nucleotide is part of comes up.
    if len(set(checklist)) == 1:
        return False,Two_Chains_Break #If length of the checklist is 1 that means there is only 1 chain, returning a false value for the Chainsbool prompting the algorithm to look for the break based on nucleotide positions 
    else:                             #Instead of taking the Chain Break based on the chains.
        return True,Two_Chains_Break

def Find_break(listed_instance):
    Pos_list=[]
    for nucleotide in listed_instance:
        split_n= nucleotide.split("|")
        Position=split_n[4]
        Pos_list.append(Position)
    for i in range(len(Pos_list)-1):
        Break = int(Pos_list[i+1]) - int(Pos_list[i])
        if Break != 1:
            Break_Index = i+1
        else:pass
    return Break_Index

def Building_sequences_ib_il(First,Second): #This Function takes a Internal Loop Instance (at this point two separated lists of the Core Nucleotides of the Motif) and extracts the bases. Also based on the length of the two parts
    #the function decides whether there is a bulge loop or and internal loop. Returns the base sequences in either case as well as a bool specifying which loop type the sequence belongs to.
    if len(First) == 0 or len(Second) == 0:
        Bulge=True
    else:
        Bulge=False
    Alpha=Building_sequences_base(First)
    Beta=Building_sequences_base(Second)
    if Bulge == True:
        sequence=Alpha+Beta
    else:
        sequence=Alpha+"$"+Beta
    return sequence,Bulge

def Building_sequences_base(loop): #Base function that takes a nucleotide list and extracts the base sequence and returns it as a string. This is used in its base form for hairpins as these are already useable like that,
    sequence=[]#In the case of Internal Loops this function is used to build the two separate Strings that are then concatenated through a $ in the function building_sequences_il_bl.
    for nucleotide in loop:
        split=nucleotide.split("|")
        sequence.append(split[3])
    sequence="".join(sequence)
    return sequence

def Output_prep(Sequence_list,id_abbreviation_dict): #Takes in a list of tuples (sequence,Loop_ID) and makes it into a set of the sequences with their single latter abbreviation with the following format: sequence+Abbreviation (e.g. GUGA+4)
    seqs=[]
    revseqs=[]
    for tuple in Sequence_list:
        String=tuple[0]+"+"+str(id_abbreviation_dict[tuple[1]])
        Rev=tuple[0][::-1]
        Rev_String=Rev+"+"+str(id_abbreviation_dict[tuple[1]])
        seqs.append(String)
        revseqs.append(Rev_String)
    seqs_set=set(seqs)
    revseqs_set=set(revseqs)
    return seqs_set,revseqs_set

if __name__ == "__main__":
    Reverse=sys.argv[1] #decides motif reverse mode, 1 only forward motifs are written, 2 only backwards motifs are written 3, both forward and reverse, 4 = all at once
    hairpin_path,internal_path,loops_path=get_filepaths() #finds the filepaths to /Data/hl_ , /Data/il_ , /Data/Loop_List
    curated_loops_dict=relevant_loops(loops_path) #Dict connecting all relevant IDs to their single letter abbreviations
    idmotif_dict=get_loop_ids(hairpin_path,internal_path) #Dict connceting all IDs to all their its motif instances
    idapi_dict=API_prep(idmotif_dict,curated_loops_dict) #Dict connecting all relevant IDs to their api calls
    idr_dict=API(idapi_dict) #Maps IDs to a List of all their API responses
    idn_dict=API_postprocess_part1(idr_dict) #Half the postprocessing done, hairpin ids now only have the motif nucleotides. Still have to get rid of the second closing basepair for internal loops in part2
    HL,IL,BL=API_postprocess_part2(idn_dict)#Finishes the postprocessing, takes out the second closing base pair for Internals and Bulge loops and sorts all sequences into their respective category (hl,il,bl) returns lists of sequence strings
    Writing_sequences(HL,IL,BL,curated_loops_dict,Reverse) #Finalizes the sequences by adding the single letter abbreviations to the strings and writes the sequences to the dedicated files in ./Loops/