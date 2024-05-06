#A simplified and class based approach to motif collection, hopefully making it less of a mess
from collections import defaultdict
import json
import os
import re
import requests
from time import sleep
from Bio import AlignIO
from io import StringIO

class Motif:
    def __init__(self,motif_json:dict):
        self.name = motif_json["motif_name"]
        self.abbreviation = motif_json["abbreviation"]
        self.instances=motif_json["instances"]
        self.loop_type=motif_json["loop_type"]
        self.rfam_ids=motif_json["rfam_id"]
        self.sequence_dict= defaultdict(list,{ key:[] for key in ["bgsu_sequences","bgsu_reverse","rfam_sequences","rfam_reverse"]})
  
    def reverse_sequences(self,seqlist:list[str]):
        reverse=[x.split(",")[0][::-1]+",{y}".format(y=x.split(",")[1]) for x in seqlist]
        return reverse
  
    def get_bgsu_sequences(self):
        for instance in self.instances:
            self.sequence_dict["bgsu_sequences"].extend(instance.get_sequences())
        
    def remove_sequence(self,sequence):
        try:self.sequence_dict["bgsu_sequences"].remove(sequence)
        except:pass
        try:self.sequence_dict["bgsu_reverse"].remove(sequence[::-1])
        except:pass
        try:self.sequence_dict["rfam_sequences"].remove(sequence)
        except:pass
        try:self.sequence_dict["rfam_reverse"].remove(sequence[::-1])
        except:pass   
          
class Hairpin(Motif):

    def __init__(self, motif_json,bgsu_json):
        super().__init__(motif_json)
        try:self.rfam_lower_bound=int(motif_json["rfam_lower_bound"])
        except:pass
        try:self.rfam_upper_bound=int(motif_json["rfam_upper_bound"])
        except:pass
        try:self.get_instances(bgsu_json)
        except:
            print("The motif {mot} is not in bgsu_json, collection will continue without this part.".format(mot=self.name))
    
    def get_instances(self,bgsu:list):
        alignments=[]
        for i in range(len(bgsu)):
            ID=re.split("[.]",bgsu[i]["motif_id"])[0]
            if ID in self.instances:
                alignments.append(Instance("hairpin",ID,bgsu[i]["alignment"],int(bgsu[i]["num_nucleotides"]),self.abbreviation))
        self.instances = alignments
    
    def get_rfam_sequences(self):
        self.rfam_api_calls  = self.make_rfam_api_calls()
        self.rmfam_alignment = self.get_rfam_alignments()
        self.sequence_dict["rfam_sequences"] += self.extract_rmfam_sequences()

    def make_rfam_api_calls(self):
        calls=[]
        for entry in self.rfam_ids:
            call="https://rfam.org/motif/{mot}/alignment?acc={mot}&format=stockholm&download=0".format(mot=entry)
            calls.append(call)
        return calls

    def get_rfam_alignments(self):
        for call in self.rfam_api_calls:
            answer=requests.get(call)
            if answer.status_code == 200:
                self.rfam_api_calls.remove(call)
                decoded=answer.content.decode()
                Alignment = list(AlignIO.parse(StringIO(decoded), format="stockholm"))[0]
                return Alignment #since all the hairpins only have one single RMFAM ID, this works. If
                                 #I ever add this to internal loops I'll have to return a list of alignments
            else:
                raise ConnectionError("Error during API request for {mot}, request return:{code}".format(mot=self.name,code=answer.status_code))
    
    def extract_rmfam_sequences(self):
        sequences=[]
        for entry in self.rmfam_alignment._records: #type:ignore
            seq=[x.upper() for x in list(entry.seq[self.rfam_lower_bound:self.rfam_upper_bound]) if x!= "-"]
            if "N" not in seq:
                sequences.append(''.join(seq)+",{abb}".format(abb=self.abbreviation))
            else:pass
        return list(set(sequences))
       
class Internal(Motif):

    def __init__(self, motif_json,bgsu_json):
        super().__init__(motif_json)
        self.get_instances(bgsu_json)

    def get_instances(self,bgsu:list):
        alignments=[]
        for i in range(len(bgsu)):
            ID=re.split("[.]",bgsu[i]["motif_id"])[0]
            if ID in self.instances:
                alignments.append(Instance("internal",ID,bgsu[i]["alignment"],int(bgsu[i]["num_nucleotides"]),self.abbreviation,int(bgsu[i]["chainbreak"])))
        self.instances=alignments

    def get_rfam_sequences(self): #Since I went through the effort of writing down every rfam internal sequence in that file, this is where they come from.
        rfam_internals_path=os.path.dirname(os.path.realpath(__file__))+ "/data/rfam_internals_fw.csv"
        with open(rfam_internals_path,"r")  as file:
            rows=file.readlines()
        split_rows=[x.split(",") for x in rows] #idk why but I have to do the splitting in a separate step, creating an extra list but the file aint big so should be no problem
        self.sequence_dict["rfam_sequences"].extend(x[0]+",{abb}".format(abb=self.abbreviation) for x in split_rows if x[1].strip() == self.abbreviation)

    def filter_sequences(self):
        keys = list(self.sequence_dict.keys())
        for key in keys:
            self.sequence_dict["i"+key] = [x for x in self.sequence_dict[key] if "$" in x]
            self.sequence_dict["b"+key] = [x for x in self.sequence_dict[key] if "$" not in x]
        for key in keys:
            del self.sequence_dict[key]

class Instance:
    def __init__(self,looptype,id:str,alignments:dict,len:int,abb:str,chainbreak:int=0):
        self.loop_type      = looptype
        self.id             = id
        self.alignments     = alignments
        self.length         = len
        self.chainbreak     = chainbreak
        self.api_call       = "http://rna.bgsu.edu/correspondence/pairwise_interactions_single?selection_type=loop_id&selection=" #full api call that every instance needs
        self.abbreviation   = abb
 
    def get_sequences(self)-> list: #Extra function, could easily put this into __init__ but to make the algorithm more readable by calling "get sequences" on all instances
       return self.get_sequences_json() + self.get_sequences_api(self.api_requests())

    def api_requests(self) ->list[list[str]]:
        api_requests  = []
        api_responses = []
        for loop in self.alignments.keys():
            request = self.api_call+loop
            api_requests.append(request)
        while len(api_requests):
            for call in api_requests:
                response=requests.get(call)
                if response.status_code == 200:
                    api_requests.remove(call)
                    decoded=response.content.decode()
                    split=decoded.split()
                    api_responses.append(split)
                else:
                    sleep(1)
        return api_responses

    def get_sequences_json(self) -> list: #returns list of sequences of all Loops in this instance, taken from the provided .json
        sequences=[]
        if self.loop_type == "hairpin":
            for loop in self.alignments.values():
                sequence = "".join( [self.get_nucleotide_element(x,3) for x in loop[1:-1]])
                if len(sequence) > 3: #Important length check for Hairpin Loops as I want to filter out 3 nucleotide non bulged UNCGs and GNRAs
                    sequences.append(sequence+",{abb}".format(abb=self.abbreviation)) #because the database is inconsistant with the UNCGs specifically, listing all their second nucleotides as bulged
        if self.loop_type == "internal":
            for loop in self.alignments.values():
                alpha = "".join( [self.get_nucleotide_element(x,3) for x in loop[1 : self.chainbreak  -1 ]] )
                omega = "".join( [self.get_nucleotide_element(x,3) for x in loop[self.chainbreak +1 : -1 ]] )
                sequence = self.FUSION(alpha,omega)
                sequences.append(sequence+",{abb}".format(abb=self.abbreviation))
        return list(set(sequences))

    def get_sequences_api(self,api:list) -> list: #returns list of sequences of all Loops in that instance, taken from the bgsu API
        sequences=[]
        r=re.compile("[|]")
        for response_list in api:
            nucleotides=[x for x in response_list if r.search(x)]
            if self.loop_type == "hairpin":
                sequence = "".join( [self.get_nucleotide_element(x,3) for x in nucleotides[1:-1]]) #it's easy for hairpins because no sequence break, what about internals tho?
                sequences.append(sequence+",{abb}".format(abb=self.abbreviation))
            if self.loop_type == "internal":
                seq_break=self.sequence_break_api(nucleotides) #seq_break is the last nucleotide of the first part of the sequence
                alpha = "".join( [self.get_nucleotide_element(x,3) for x in nucleotides [ 1 : seq_break   ]])
                omega = "".join( [self.get_nucleotide_element(x,3) for x in nucleotides [ seq_break+2 : -1]])
                sequence = self.FUSION(alpha,omega)
                sequences.append(sequence+",{abb}".format(abb=self.abbreviation))
        return list(set(sequences))

    def sequence_break_api(self,nucleotides:list) -> int:
        positions=[ self.get_nucleotide_element(x,4) for x in nucleotides ]
        for i in range(len(positions)-1):
            Diff = int(positions[i+1]) - int(positions[i])
            if abs(Diff) >= 3: #This is where bulge size is theoretically limited, if there are two internal loop parts that are
                return i
            else:pass
        raise IndexError("There is no sequence break with 6 or more nucleotides between the internal strands")

    def FUSION(self,a,b) -> str:
        if len(a) > 0 and len(b) > 0:
            sequence = a + "$" + b #type:str
        else:
            sequence = a + b #concatenate them together, this way it doesnt matter which one has length 0. Later I can check for $ in the string to decide if its bulge or internal.
        return sequence

    def get_nucleotide_element(self,nucleotide:str,number:int) ->str:
        split=nucleotide.split("|")        
        element=split[number]
        return element

def load_local_json(file_name:str) -> list:
    local_file = os.path.dirname(os.path.realpath(__file__)) + "/data/" + file_name
    with open(local_file) as json_file:
        motifs_json=json.load(json_file)
    return motifs_json
    
def load_bgsu_json(call:str,alternative:str) ->list: #Tries to download the latest json from bgsu 5 times, if all 5 attempts fail it takes a locally provided 3.81 version from the data folder.
    i = 0
    while i < 6:
        response= requests.get(call)
        if response.status_code == 200:
            return json.loads(response.content.decode())
        else:
            i += 1
    print("Could not establish connection to BGSU servers, resorting to hl_3.81.json and il_3.81.json backup files in /data/ folder.")
    local = load_local_json(alternative)
    return local

def load_jsons() -> list[Hairpin | Internal]: #The load_bgsu function has a backup hl_3.81/il_3.81 in case the most up to date versions can not be fetched from the bgsu servers.
    hl=load_bgsu_json("http://rna.bgsu.edu/rna3dhub/motifs/release/hl/current/json","hl_3.81.json")
    il=load_bgsu_json("http://rna.bgsu.edu/rna3dhub/motifs/release/il/current/json","il_3.81.json")
    motif_json=load_local_json("motifs.json")
    motifs=[] #type:list
    for motif in motif_json:
        if motif["loop_type"] == "hairpin":
            class_motif = Hairpin(motif,hl)
            motifs.append(class_motif)
        else:pass
        if motif["loop_type"] == "internal":
           class_motif = Internal(motif,il)
           motifs.append(class_motif)
        else:pass
    return motifs

def dupe_check(looptype:str,motif_list:list[Hairpin|Internal])-> bool: #quick function that prints out doubled sequences (appearing for two different motifs).
    seen={}
    check_bool=False
    for motif in motif_list:  #As of 25.04.24 BGSU Version 3.81 there are two: UUCAA (GNRA/T-Loop) and GUGA (GNRA/UNCG)
        if motif.loop_type==looptype:
            for sequence in list(set(flatten([list(set(motif.sequence_dict["bgsu_sequences"])),list(set(motif.sequence_dict["rfam_sequences"]))]))):
                if sequence[:-2] not in seen.keys():
                    seen[sequence[:-2]]=motif.abbreviation
                else:
                    print("Duplicated sequence found: {seq}. Found sequence in {a} and {b}. If you proceed with the created header only the first instance will be recognized. Consider removing one of the two with the motif.remove_sequence function".format(seq=sequence[:-2],a=seen[sequence[:-2]],b=motif.abbreviation))
    return check_bool
           
def create_hexdumbs(motif_list:list[Hairpin|Internal]):
    keys=["bgsu_fw","bgsu_rv","bgsu_both","rfam_fw","rfam_rv","rfam_both","both_fw","both_rv","both_both"] #type:list[str]
    hsequence_dict = defaultdict(list,{ key:[] for key in keys })
    isequence_dict = defaultdict(list,{ key:[] for key in keys })
    bsequence_dict = defaultdict(list,{ key:[] for key in keys })
    for mot in motif_list:
        if isinstance(mot,Hairpin):
            sort_sequences(mot,hsequence_dict)
        if isinstance(mot,Internal):
            mot.filter_sequences()
            sort_sequences(mot,isequence_dict,"i")
            sort_sequences(mot,bsequence_dict,"b")
    with open("mot_header.hh","w+") as file:
        for key in keys:
            file.write(sequences2header(hsequence_dict[key],"h"+key))
            file.write(sequences2header(isequence_dict[key],"i"+key))
            file.write(sequences2header(bsequence_dict[key],"b"+key))

def sort_sequences(motif:Hairpin|Internal,seq_dict:dict[str,list[str]],looptype:str=""): #modifies the dictionary inplace, so no return needed
    seq_dict["bgsu_fw"].extend(list(set(motif.sequence_dict[looptype+"bgsu_sequences"])))
    seq_dict["bgsu_rv"].extend(list(set(motif.sequence_dict[looptype+"bgsu_reverse"])))
    seq_dict["bgsu_both"].extend(list(set(flatten([motif.sequence_dict[looptype+"bgsu_sequences"],motif.sequence_dict[looptype+"bgsu_reverse"]]))))
    seq_dict["rfam_fw"].extend(list(set(motif.sequence_dict[looptype+"rfam_sequences"])))
    seq_dict["rfam_rv"].extend(list(set(motif.sequence_dict[looptype+"rfam_reverse"])))
    seq_dict["rfam_both"].extend(list(set(flatten([motif.sequence_dict[looptype+"rfam_sequences"],motif.sequence_dict[looptype+"rfam_reverse"]]))))
    seq_dict["both_fw"].extend(list(set(flatten([motif.sequence_dict[looptype+"bgsu_sequences"],motif.sequence_dict[looptype+"rfam_sequences"]]))))
    seq_dict["both_rv"].extend(list(set(flatten([motif.sequence_dict[looptype+"bgsu_reverse"],motif.sequence_dict[looptype+"rfam_reverse"]]))))
    seq_dict["both_both"].extend(list(set(flatten([flatten([motif.sequence_dict[looptype+"bgsu_sequences"],motif.sequence_dict[looptype+"rfam_sequences"]]),flatten([motif.sequence_dict[looptype+"bgsu_reverse"],motif.sequence_dict[looptype+"rfam_reverse"]])]))))

def flatten(xss:list[list[str]]) -> list:
    return [x for xs in xss for x in xs]

def sequences2header(seq_set:list,name:str)->str:
    joined_seq_set="\n".join(seq_set)
    out=[]
    out.append('static unsigned char {var_name}[] = {{'.format(var_name=name))
    data =[ joined_seq_set[ i : i+12 ] for i in range(0, len(joined_seq_set), 12) ]
    for i,x in enumerate(data):
        line=', '.join((["0x{val:02x}".format(val=ord(c)) for c in x]))
        out.append('  {lined}{comma}'.format(lined=line,comma=',' if i <len(data)-1 else ''))
    out.append('};')
    out.append('static unsigned int {var_name}_len = {data_len};\n'.format(var_name=name,data_len=len(joined_seq_set)))
    return '\n'.join(out)

if __name__ == "__main__":
    motifs=load_jsons()
    for motif in motifs:
        if len(motif.instances) > 0:
            motif.get_bgsu_sequences()
            motif.sequence_dict["bgsu_reverse"]=motif.reverse_sequences(motif.sequence_dict["bgsu_sequences"])
        if len(motif.rfam_ids) > 0: #Just a check for whether the motif is also present in Rfam, indicated by whether there are any rfam ids given in the motif.json file
            motif.get_rfam_sequences()
            motif.sequence_dict["rfam_reverse"]=motif.reverse_sequences(motif.sequence_dict["rfam_sequences"])
        #removes the two currently dubed sequences, this needs to be done manually as the program has no idea what motif the sequence should belong to!
        if motif.name == "UNCG":
            motif.remove_sequence("GUGA,U")
        if motif.name == "GNRA":
                motif.remove_sequence("UUCAA,G")
    #Here all the sequences need to be unambiguous, so potetially doubled sequences need to be removed during the for loop by checking the motif.name and then remove the entry with the input type being a string formatted: sequence,abbreviation.
    Check1=dupe_check("hairpin",motifs)
    Check2=dupe_check("internal",motifs)
    if not Check1 and not Check2:
        create_hexdumbs(motifs)
    else:
        raise NameError("Hexdumbs were not created because sequences are not unambiguous. Please delete sequences with remove_sequence function")