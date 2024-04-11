##     script for analyse of exonerate outputs looking for NLR genes
##
##   by Lucas Bonometti

###############################################################################
###############################################################################
########
########                  libraries
########
###############################################################################
###############################################################################

import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
from itertools import combinations
from timeit import default_timer as timer
from datetime import timedelta

###############################################################################
###############################################################################
########
########                  the parameters of the file
########
###############################################################################
###############################################################################

parser = argparse.ArgumentParser()

parser.add_argument("-n","--fasta_allnlrs",type=str,required=False,
                    help="the fasta containing all the NLRs, as nucleotide sequences") 

parser.add_argument("-d","--dir_ogs",type=str,required=False,
                    help="directory containing the nucleotide fastas of all the NLR orthogroups,"+
                    " with both the NLRs and the non-NLRs genes of those orthogroups")

parser.add_argument("-g","--gff_file",type=str,required=False,
                    help="the gff file associated to the assembled genome")

parser.add_argument("-a","--assembly_file",type=str,required=True,
                    help="the assembly file of the genome")

parser.add_argument("-e","--exo_all",type=str,required=True,
                    help="the file containing the exonerate output for the run of all the NLRs"+
                    " on the analyzed genome, with gff included")

parser.add_argument("-b","--exo_nb",type=str,required=True,
                    help="the file containing the exonerate output for the run of all the NBs"+
                    " on the analyzed genome, without gff included")

parser.add_argument("-t","--threshold_nlr",type=int,required=False,default=100,
                    help="threshold for accepting an exonerate hit of a NLR as a valide hit "+
                    "(default = 100)")

parser.add_argument("-r","--threshold_nb",type=int,required=False,default=100,
                    help="threshold for accepting an exonerate hit of a NB as a valide hit "+
                    "(default = 100)")

parser.add_argument("-m","--minsize_nlr",type=int,required=False,default=1000,
                    help="mininal size, in bp, for accepting a gene or a hit as a valide NLR"+
                    "(default = 1000)")

parser.add_argument("-c","--perc_shortened_cterm",type=int,required=False,default=10,
                    help="proportion, in percent, of a query gene that can be shorten in "+
                    "the c-terminal extremity of the target hit (default = 10)")

parser.add_argument("-f","--fasta_nlrs_nt_out",type=str,required=False,
                    help="output fasta of all the NLRs finally identified in the given genome, in ucleotides")

parser.add_argument("-z","--fasta_nlrs_aa_out",type=str,required=False,
                    help="output fasta of all the NLRs finally identified in the given genome, in amino acids")

parser.add_argument("-p","--pseudo_out",type=str,required=False,
                    help="output table of all the pseudo-NLRs identified in the given genome")

parser.add_argument("-o","--output",type=str,required=False,
                    help="output file containing all the NLRs and pseudo-NLRs identified,"+
                    " and their respective gffs")

args = parser.parse_args()

###############################################################################
###############################################################################
########
########                  the functions
########
###############################################################################
###############################################################################


###############################################################################
######## function 1 : 
######## takes a gff as entry and returns
######## the sum of the lengths of its exons
###############################################################################

def longueur(gff):
    l=0
    for x in gff.keys():
        l=l+1+max(gff[x])-min(gff[x])
    return l

###############################################################################
######## function 2 : 
######## takes a gff and its frame ("+" or "-") as entry
######## and returns the sorted gff
###############################################################################

def sortgff(gff,f):
    inef=0
    newgff={}
    
# we sort in two different ways depending on wether the gene in reverse-complement or not
    if f=="+":
        revsens=False
    elif f=="-":
        revsens=True
        
# we sort the gff based on the min value of each exon (min, max)
# we then suppose the exons are non-overlapping
    for kmin in sorted({min(kj) for kj in gff.values()},reverse=revsens):
        for ji in gff.values():
            if kmin==min(ji):
                newgff[inef]={kmin,max(ji)}
                inef+=1
                break
    return newgff

###############################################################################
######## function 3 : 
######## takes a gff, its frame and its associated scaffold (as a seqio record) as entry
######## returns true if there is no stop codon but at the last position
######## returns false otherwise
###############################################################################

def gffisok(gff,frame,scarecord):
    ntsq=""
    gf=sortgff(gff,frame)
    for i in range(len(gf)):
        if len(gf[i])>1:
            sq=scarecord.seq[min(gf[i])-1:max(gf[i])]
        else:
            sq=Seq(scarecord.seq[min(gf[i])-1])
        if frame=="-":
            sq=sq.reverse_complement()
        ntsq=ntsq+str(sq)
    protsq=str(Seq(ntsq).translate())
    if "*" in protsq[:-1]:
        return False
    else:
        return True

###############################################################################
######## function 4 : 
######## takes two gffs (gff1 and gff2), their common frame ("+" or "-")
######## and the length of their scaffold as entry
######## and returns a new hybrid gff
######## gff1 is priorized over gff2
###############################################################################

def mergegff(gff1,gff2,frame,scarecord):
    
# mingf1 and maxgf1 are respectively the minimal and maximal positions of gff1
    mingf1= min({min(kj) for kj in gff1.values()})
    maxgf1= max({max(kj) for kj in gff1.values()})
    
# mingf2 and maxgf2 are respectively the minimal and maximal positions of gff2
    mingf2= min({min(kj) for kj in gff2.values()})
    maxgf2= max({max(kj) for kj in gff2.values()})
    
# (1) first situation
# as gff1 is priorized over gff2, if gff2 is included in the limits of gff1
# then gff1 is returned
    if mingf1<=mingf2 and maxgf1>=maxgf2:
        
        return gff1
    
# if the first situation is not reached
# then two cases are considered, one for each possible farme ("+" or "-")
# only the code for the "+" frame will be commented, the code for the "-" frame mirroring it
    elif frame=="+":
        
# gf1 is the sorted gff1 and gf2 is the sorted gff2
        gf1=sortgff(gff1,frame)
        gf2=sortgff(gff2,frame)
        
# frame1 is the frame (0, 1 or 2) of each element of gf1 and frame2 of gf2
# ex. frame1[0] is the frame of gf1[0] and frame1[3] of gf1[3]
        frame1={}
        frame2={}
        
# multi3 is True if the sequence needs to be a multiple of 3, False otherwise
# multi3 is False only if the sequence end is at a scaffold extremity
        multi3=True
        
# f is the current frame (0, 1 or 2)
        f=0

# this part implements frame1 and frame2
        for i1 in sorted(gf1.keys()):
            frame1[i1]=f
            f=(f+(max(gf1[i1])-min(gf1[i1])+1)%3)%3
        f=0
        for i2 in sorted(gf2.keys()):
            frame2[i2]=f
            f=(f+(max(gf2[i2])-min(gf2[i2])+1)%3)%3

# here we check if the gff1 ends at an extremity of the scaffold, and the same for gff2
# if it is the case, we considere the gene as possibly truncated,
# and then  we allow it to not be a multiple of 3
        if max({max(kj) for kj in gf1.values()})==len(scarecord.seq):
            multi3=False
        if max({max(kj) for kj in gf2.values()})==len(scarecord.seq):
            multi3=False
        
# toaddbefore and toaddafter are dictionnaries of exons of gf2 to add to gf1,
# respectively before the begining and after the end of gf1
# frameadbefore and frameadafter are the dictionnary of the frames of the related pieces
# it is the same for toaddcutbefore and toaddcutafter, and their farmes frameadcutbefore and frameadcutafter,
# but there the exons from gf2 are truncated because overlapping with the first or last exon of gf1
        toaddbefore={}
        frameadbefore={}
        toaddcutbefore={}
        frameadcutbefore={}
        toaddafter={}
        frameadafter={}
        toaddcutafter={}
        frameadcutafter={}
        
# iad is the rank in toadd and in framead, with :
# a=after
# b=before
# c=cut
# ex. if iadb=2 it means that frameadbefore (and then frameadbefore too) contains 2 elements (2 exons)
        iadb=0
        iadcb=0
        iada=0
        iadca=0
        
# this for loop is made to take the exons from gf2 and add them, or not, in the toadd dicos
# the frame is then estimated from the frame of the whole exon of gf2,
# except when it is a cut-after exon
        for i2 in gf2.keys():
            
# there if the given exon from gf2 is before gf1
            if min(gf2[i2])<mingf1:
                
# and there if the given exon from gf2 overlaps with gf1
                if max(gf2[i2])>=mingf1:
                    toaddcutbefore[iadcb]={min(gf2[i2]),mingf1-1}
                    frameadcutbefore[iadcb]=frame2[i2]
                    iadcb+=1

# or there if it does not
                else:
                    toaddbefore[iadb]=gf2[i2]
                    frameadbefore[iadb]=frame2[i2]
                    iadb+=1

# there is if the given exon from gf2 is after gf1
            if max(gf2[i2])>maxgf1:
                
# and there if the given exon from gf2 overlaps with gf1
# the frame of the truncated exon is then recalculated based on the frame of the whole exon of gf2
                if min(gf2[i2])<=maxgf1:
                    toaddcutafter[iadca]={maxgf1+1,max(gf2[i2])}
                    frameadcutafter[iadca]=(maxgf1-min(gf2[i2])+1+frame1[i1])%3
                    iadca+=1

# or there if it does not
                else:
                    toaddafter[iada]=gf2[i2]
                    frameadafter[iada]=frame2[i2]
                    iada+=1
# liset is the list of the toadd dicos and their frames dicos
# if the toadd dicos are not empty, their rank (0, 1, 2 or 3) are added to comblist
        liset=[[toaddbefore,frameadbefore],
               [toaddcutbefore,frameadcutbefore],
               [toaddcutafter,frameadcutafter],
               [toaddafter,frameadafter]]
        comblist=[]
        for icomb in range(4):
            if len(liset[icomb][0])!=0:
                comblist.append(icomb)
        
# comb is a list of the combination of all the non-empty keys of comblist
# from all the keys (ex. [0,1,2,3] to only one (ex. [2])
        comb=[]
        for nbad in range(4,0,-1):
            comb=comb+list(combinations(comblist,nbad))
        
# here we will add the pieces of the toadd dicos to gf1 and see if two conditions are respected :
# 1. newgff needs to be a multiple of 3 (or not if multi3 is False)
# 2. all the frames of the exons must be respected
# to do so, we test the different combinations of comb, beginning with the ones with most exons
        for c in comb:
            
# newgff is gf1 with the new pieces from gf2, and newframe its related frames
            newgff=gf1
            newframe=frame1
            
# this for loop is used to add the new exons to newgff
            for ic in c:
                for isadded in range(len(liset[ic][0])):
                    newgff[len(newgff)]=liset[ic][0][isadded]
                    newframe[len(newgff)-1]=liset[ic][1][isadded]
            
# there we check the first condition (= is a multiple of 3)
            if multi3==True:
                m3=longueur(newgff)%3==0
            elif multi3==False:
                m3=True
            
# there we check the second condition (= all frames are correct)
# to do so, first we sort newgff as newgf,
# and we create the related correct (because based on newframe) frames dico (fr)
            goodframe=True
            newgf=sortgff(newgff,frame)
            fr={}
            for ine in newgf.keys():
                for ike in newframe.keys():
                    if newgf[ine]==newgff[ike]:
                        fr[ine]=newframe[ike]
            
# then we test if following newgf we have the same frame f as expected in fr
            f=0
            for ine in sorted(newgf.keys()):
                if f!=fr[ine]  : goodframe=False
                f=(f+(max(newgf[ine])-min(newgf[ine])+1)%3)%3
            
# (2) if both conditions are respected, we stop here and we return newgf
# (3) if at least one condition is never respected, we return gff1
            if m3==True and goodframe==True:
                break
        if m3==False or goodframe==False:
            return gff1
        else:
            if gffisok(newgf,frame,scarecord):
                return newgf
            else:
                return gff1

# this part is the same as for frame "+" but with frame "-"
# consequently, some tiny differences are present in the code
# but the base is identical
    elif frame=="-":
        
        gf1=sortgff(gff1,frame)
        gf2=sortgff(gff2,frame)
        
        frame1={}
        frame2={}
        multi3=True
        f=0
        for i1 in sorted(gf1.keys()):
            frame1[i1]=f
            f=(f+(max(gf1[i1])-min(gf1[i1])+1)%3)%3

# difference here : the limit of the scaffold is position "1", not its length
        if min({min(kj) for kj in gf1.values()})==1:
            multi3=False
        f=0
        for i2 in sorted(gf2.keys()):
            frame2[i2]=f
            f=(f+(max(gf2[i2])-min(gf2[i2])+1)%3)%3
        if min({min(kj) for kj in gf2.values()})==1:
            multi3=False
        
        toaddbefore={}
        frameadbefore={}
        toaddcutbefore={}
        frameadcutbefore={}
        
        toaddafter={}
        frameadafter={}
        toaddcutafter={}
        frameadcutafter={}
        
        iadb=0
        iadcb=0
        iada=0
        iadca=0
        for i2 in gf2.keys():
            if min(gf2[i2])<mingf1:
                if max(gf2[i2])>=mingf1:
                    toaddcutbefore[iadcb]={min(gf2[i2]),mingf1-1}

# difference here : the calculation of the frame is different 
                    frameadcutbefore[iadcb]=(max(gf2[i2])-mingf1+1+frame1[i1])%3
                    iadcb+=1
                else:
                    toaddbefore[iadb]=gf2[i2]
                    frameadbefore[iadb]=frame2[i2]
                    iadb+=1
            if max(gf2[i2])>maxgf1:
                if min(gf2[i2])<=maxgf1:
                    toaddcutafter[iadca]={maxgf1+1,max(gf2[i2])}
                    frameadcutafter[iadca]=frame2[i2]
                    iadca+=1
                else:
                    toaddafter[iada]=gf2[i2]
                    frameadafter[iada]=frame2[i2]
                    iada+=1
        
        liset=[[toaddbefore,frameadbefore],
               [toaddcutbefore,frameadcutbefore],
               [toaddcutafter,frameadcutafter],
               [toaddafter,frameadafter]]
        
        comblist=[]
        for icomb in range(4):
            if len(liset[icomb][0])!=0:
                comblist.append(icomb)
        
        comb=[]
        for nbad in range(4,0,-1):
            comb=comb+list(combinations(comblist,nbad))
        
        for c in comb:
            newgff=gf1
            newframe=frame1
            for ic in c:
                for isadded in range(len(liset[ic][0])):
                    newgff[len(newgff)]=liset[ic][0][isadded]
                    newframe[len(newgff)-1]=liset[ic][1][isadded]
            if multi3==True:
                m3=longueur(newgff)%3==0
            elif multi3==False:
                m3=True
            goodframe=True
            newgf=sortgff(newgff,frame)
            fr={}
            for ine in newgf.keys():
                for ike in newframe.keys():
                    if newgf[ine]==newgff[ike]:
                        fr[ine]=newframe[ike]
            f=0
            for ine in sorted(newgf.keys()):
                if f!=fr[ine]  : goodframe=False
                f=(f+(max(newgf[ine])-min(newgf[ine])+1)%3)%3
                
            if m3==True and goodframe==True:
                break
        if m3==False or goodframe==False:
            return gff1
        else:
            if gffisok(newgf,frame,scarecord):
                return newgf
            else:
                return gff1

# this is set to print an error in the case the frames are not correct
    else:
        print("problem of frame : neither '-' or '+'")
        return gff1

###############################################################################
###############################################################################
########
########              intiation of statistical parameters
########
###############################################################################
###############################################################################

#
sacers_ind=0

#
sacer_notkept=set()

#
sacer_kept=0

#
sacer_nlrnonb=0

#
nlr_ind=0

#
nlr_longer=0

#nlr_notouch is the number of NLR genes that are finally not affected by the correction
nlr_notouch=0

#
nlr_nonb=0

#
newnlr_nosacer=0

#
newnlr_unknown=0

#
unknown_nonb=0

#
nlrlist_ind=set()

#
nlrlist_prolong=set()

# nlrlist_notouch is dico of all the NLR genes that are not affected by the correction,
# and their scaffold
nlrlist_notouch={}

#
nlrlist_pseudo=set()

#
sacerlist_kept=set()

#
newgenelist=set()

# t0 is the beginning running time, as reference for later on
t0=timer()

# ind is the code for the individual (ex. ind=neucra-01)

ind=args.exo_all.split("/")[-1].split("_")[2].split(".")[0]

###############################################################################
###############################################################################
########
########                 the core of the script
########
###############################################################################
###############################################################################

###############################################################################
######## 1. introduction
###############################################################################

# readable, here and latter, is a readable output, explaining how the script has run
with open("./readable_"+ind+"_correc_exo.txt",'w') as readable:
    
    readable.write("#START\n")
    readable.write("#\n")
    
# assembly is a dico of the scaffolds and their length in number of nucleotides
    assembly={}
    for record in SeqIO.parse(args.assembly_file,"fasta"):
        assembly[str(record.id)]=len(record.seq)
    
# nlrs is a dico of all the old nlrs and their respective length in nucleotides
# nlr_new is a set of the nlrs with a name finishing by "_new", because manually corrected
# the "new" nlrs are latter directly considered as "bad" so that their old gff is no more considered
    nlrs={}
    nlr_new=set()
    if args.fasta_allnlrs:
        for record in SeqIO.parse(args.fasta_allnlrs,"fasta"):
            nlrs["_".join(str(record.id).split("_")[:2])]=len(record.seq)*3
            if str(record.id).split("_")[-1]=="new":
                nlr_new.add("_".join(str(record.id).split("_")[:2]))
            if str(record.id).split("_")[0]==ind:
                nlr_ind+=1
                nlrlist_ind.add("_".join(str(record.id).split("_")[:2]))

    readable.write("#exonerate analysis realized based on a set of "+str(len(nlrs))+" NLR genes\n")
    
# sacers is a dico of all the non-nlr genes of the NLR-containing OGs and their respective lengths
    sacers={}
    if args.dir_ogs:
        for file in os.listdir(args.dir_ogs):
            for record in SeqIO.parse(args.dir_ogs+file,"fasta"):
                if "_".join(str(record.id).split("_")[:2]) not in nlrs:
                    sacers["_".join(record.id.split("_")[:2])]=len(record.seq)
    
    readable.write("#"+str(len(sacers))+" non-NLR genes on the orthogroups of NLRs ('sacer' genes)\n")
    readable.write("#\n")
    readable.write("#analyse of the genomes of "+ind+"\n")
    readable.write("#start of the analyze of the NB exonerate output\n")
    readable.write("#time is "+str(timedelta(seconds=(timer()-t0)))+"\n")
    readable.write("#\n")

###############################################################################
######## 2. checking the exonerate output for "NB" domains
###############################################################################

# "nb" represents all the hits for a NB
# with information about the scaffold, start, end, score and status (good or not)
# a good NB is a sequence without any stop codon or any frame shift in it
nb={}

# here is read the exonerate output, directly, line by line
with gzip.open(args.exo_nb,'rt') as IN:
    
# we initialize the values of the hits
# scaffold is the scaffold
    scaffold=None
    
# rawscore is the exonerate score of the hit
    rawscore=None

# debut is the begining position of the hit
    debut=None
    
# fin is the ending position of the hit
    fin=None

# status is where the hit is good (no frame shift and no stop codon) or bad
    status=None

# inseq is true if we are reading the alignment, false otherwise
    inseq=False
    
    for line in IN:
        if line[0]=="#":continue
        if len(line.strip())>0:
            
# elet[0] is the class of the result given for the hit
            elet=line.strip().split(":")

# target gives the scaffold
            if elet[0]=="Target":
                scaffold=elet[1].split()[0].strip()
            
# raw score give sthe score of the exonerate hit for the NB
            elif elet[0]=="Raw score":
                rawscore=elet[1].strip()
            
# target range gives the limits of the hit on the target scaffold
# as it is the last result given before the alignment, it is used for initiating the lecture of the alignment
            elif elet[0]=="Target range":
                debut=min(int(elet[1].strip().split()[0]),int(elet[1].strip().split()[2]))
                fin=max(int(elet[1].strip().split()[0]),int(elet[1].strip().split()[2]))

# inseq is true if we read the alignment, false otherwise
                inseq=True
                
# cl is the line number as we read in the alignment
                cl=0

# status is good if there is no stop codon and no frame shift, bad otherwise
                status = "good"

# vulgar gives the information that the alignment lecture is over, so inseq in set to false
            elif elet[0]=="vulgar":
                inseq=False

# we considere only hits with score above threshold
# we then implement nb
                if int(rawscore)>=int(args.threshold_nb):
                    if scaffold not in nb:
                        nb[scaffold]={}
                        nb[scaffold][0]={"debut":debut,
                                         "fin":fin,
                                         "score":int(rawscore),
                                         "status":status}
                    else:
                        
# we compare the analyzed hit with hits already implemented in nb,
# and we keep only one hit when hits are overlapping, prefering good hits over bad hits
                        hitlist=set(nb[scaffold]).copy()
                        overlaps=False
                        for i in hitlist:
                            if nb[scaffold][i]["debut"]<fin and nb[scaffold][i]["fin"]>debut:
                                overlaps=True
                                if (status=="good" and nb[scaffold][i]["status"]=="good") or \
                                    (status=="bad" and nb[scaffold][i]["status"]=="bad"):
                                    nb[scaffold][i]["debut"] = min(debut,nb[scaffold][i]["debut"])
                                    nb[scaffold][i]["fin"] = max(fin,nb[scaffold][i]["fin"])
                                    nb[scaffold][i]["score"] = max(int(rawscore),nb[scaffold][i]["score"])
                                elif status=="good" and nb[scaffold][i]["status"]=="bad":
                                    del nb[scaffold][i]
                                    if len(nb[scaffold])>0:
                                        nb[scaffold][max(nb[scaffold].keys())+1]={"debut":debut,
                                                                                  "fin":fin,
                                                                                  "score":int(rawscore),
                                                                                  "status":status}
                                    else:
                                        nb[scaffold][0]={"debut":debut,
                                                         "fin":fin,
                                                         "score":int(rawscore),
                                                         "status":status}
                        if overlaps==False:
                            nb[scaffold][max(hitlist)+1]={"debut":debut,
                                             "fin":fin,
                                             "score":int(rawscore),
                                             "status":status}

# when we read the alignment, we only check that there is no frame shift and no stop codon
# otherwise, the hit is classified as having a bad status
            elif inseq==True:
                cl+=1
            
# cl%5=4 when we read the amino acid line of the target
# cl%5=0 when we read the nucleotide line of the target
                if cl%5==4 or cl%5==0:
                    if "*" in line or "#" in line:
                        status = "bad"

# after the lecture of the exonerate output, we merge the overlapping NB hits
# (it was done at the previous step, but some hits may have been missed)
for scaffold in nb:

# the while loop is used for the merging
# every time we merge two hits, cleared is set to false and the loop begin again
    cleared=False
    while cleared==False:
        cleared=True

# to break is false if the while loop needs to be begun again
        tobreak=False
        hitlist=set(nb[scaffold]).copy()
        for i in hitlist:
            if tobreak==True:
                break
            
# we compare every hit i to every other hit j
            for j in hitlist:
                if j==i:continue
                
# we merge overlapping hits, prefering good hits to bad hits
                if nb[scaffold][i]["debut"]<nb[scaffold][j]["fin"] and \
                    nb[scaffold][j]["debut"]<nb[scaffold][i]["fin"]:
                    if (nb[scaffold][j]["status"]=="good" and nb[scaffold][i]["status"]=="good") or \
                        (nb[scaffold][j]["status"]=="bad" and nb[scaffold][i]["status"]=="bad"):
                        nb[scaffold][i]["debut"] = min(nb[scaffold][j]["debut"],nb[scaffold][i]["debut"])
                        nb[scaffold][i]["fin"] = max(nb[scaffold][j]["fin"],nb[scaffold][i]["fin"])
                        nb[scaffold][i]["score"] = max(nb[scaffold][j]["score"],nb[scaffold][i]["score"])
                        del nb[scaffold][j]
                    elif nb[scaffold][j]["status"]=="good" and nb[scaffold][i]["status"]=="bad":
                        del nb[scaffold][i]
                    elif nb[scaffold][j]["status"]=="bad" and nb[scaffold][i]["status"]=="good":
                        del nb[scaffold][j]
                    cleared=False
                    tobreak=True
                    break

with open("./readable_"+ind+"_correc_exo.txt",'a') as readable:    
    readable.write("#NB exonerate output analyzed\n")
    readable.write("#time is "+str(timedelta(seconds=(timer()-t0)))+"\n#\n")
    readable.write("#\n#analyse of the whole-NLRs exonerate output\n")

###############################################################################
######## 3. checking the exonerate output for all the NLR genes
###############################################################################

# exhits is the dico of all the non-overlapping hits with no frame shift and no stop codon inside
exhits={}

# pseudonlr is the dico of the non-overlapping hits with at least one frame shift or one stop codon inside
# hits with a bad NB domain are directly included inside
# other are first transfered to pseudopossible
pseudonlr={}

# pseudopossible is the dico of the hits with at least one frame shift or one stop codon inside
# but not in the NB domain, so a trucated version of themselves could be included in exhits
# if not, they will be tranfered to pseudonlr
# pseudopossible is then a transitory dictionnary
pseudopossible={}

# here is read the exonerate output, directly, line by line
with gzip.open(args.exo_all,"rt") as IN:
    
# query is the name of the nlr originating the hit
    query=None
    
# scaffold is the scaffold of the hit
    scaffold=None
    
# rawscore is the exonerate score of the hit
    rawscore=None
    
# debut is the start of the hit in the scaffold
    debut=None
    
# fin is the end of teh hit in the scaffold
    fin=None
    
# qlength is the length of the query matching the hit, in nt
    qlength=None
    
# isgff is true if we are reading a gff line and false otherwise
    isgff=False
    
# sens is the frame ("+" or "-")
    sens=None
    
# gffexo is the gff from exonerate
    gffexo={}
    
# alignt is true if we are reading a line in the alignment between the query and the target,
# false otherwise
    alignt=False
    
# cl is the count of the lines in the alignment
    cl=0
    
# cdsaa is a list containing the sequences of the exons of the target gene, in the order, in amino acids
# it is based on the sequences of the alignment
    cdsaa=[""]
    
# cdsnt is the same as cdsaa, but in nucleotides
    cdsnt=[""]
    
# totreat is true if the hits overlaps with a NB domain, false otherwise
# the hit is kept only if totreat is true
    totreat=False
    
# ispseudo1 is true if the hit overlaps with a "bad" NB domain, false otherwise
# if true, the hit is included in pseudonlr
    ispseudo1=False
    
# ispseudo2 is true if the hit overlaps with a "good" NB domain but also includes at least one frame shift or
# one stop codon inside, elsewhere
# otherwise, it is false
# if true, the hit is included in pseudopossible
    ispseudo2=False
    
    for line in IN:
        
# this line means we are analyzing an new hit
        if line[:13]=="C4 Alignment:":
            totreat=True
            alignt=False
            isgff=False
            gffexo={}

# this line means we are begining to read the exonerate gff
        elif line.strip()=="# --- START OF GFF DUMP ---" and totreat == True and ispseudo1 == False:     
            isgff=True
            
# this line means that we have finished to read the gff, and then the whole hit is read
# consequently, the hit is analyzed
        elif line.strip()=="# --- END OF GFF DUMP ---" and totreat == True and ispseudo1 == False:
            isgff=False

            gffexo=sortgff(gffexo,sens)
            
# this is a test to check if we have corectly read both alignment and gff
            if len(cdsaa)!=len(gffexo.keys()) or len(cdsnt)!=len(gffexo.keys()):
                print("problem of not correspondance between gff and cds lists")
                print(cdsaa)
                print(cdsnt)
                print(gffexo)
                os.kill()

# we keep a hit only if its score is above threashold
            if int(rawscore)>=int(args.threshold_nlr):           
                
# we analyze if the gene contains stop codons (*) or frameshilfts (#)
# ispseudo2 is then set to true
# there is a tolerance for the last codon of the last exon, which can be a stop codon
# there is a tolerance in the c-terminal part of the gene, depending on parameter -c
# but only on the last exon
                ispseudo2=False
                rangexon=0

# shortinc is the number of nucleotides the hit can be reduced at most in its c-terminal extremity
# of the last exon
# it is calculated based on the query length, the hit length, and the value of the parameter -c
                shortinc=max(0,-(nlrs[query]-qlength-int(nlrs[query]*args.perc_shortened_cterm/100)))
                for exonaa in cdsaa:
                    rangexon+=1

# if we are working on the last exon, we can reduce its length if there is a stop codon or a frame shift
# after the limit set by shortinc
                    if rangexon==len(cdsaa):
                        if len(exonaa)>shortinc:
                            if "*" in exonaa[:-(3+shortinc)]:
                                ispseudo2=True
                            if shortinc!=0:
                                if "#" in exonaa[:-shortinc]:
                                    ispseudo2=True
                            elif "#" in exonaa:
                                ispseudo2=True
                        if ispseudo2==False:
                            if "*" in exonaa[:-3] or "#" in exonaa:
                                lastpos=min(exonaa.find("*"),exonaa.find("#"))-1
                                toshorten=len(exonaa)-lastpos+1
                                if sens=="+":
                                    if max(gffexo[max(gffexo.keys())])-toshorten < \
                                        min(gffexo[max(gffexo.keys())]):
                                            gffexo.pop(max(gffexo.keys()))
                                    else:
                                        gffexo[max(gffexo.keys())]={min(gffexo[max(gffexo.keys())]),
                                                                    max(gffexo[max(gffexo.keys())])-toshorten}
                                    gffexo=sortgff(gffexo,sens)
                                elif sens=="-":
                                    if min(gffexo[max(gffexo.keys())])+toshorten > \
                                        max(gffexo[max(gffexo.keys())]):
                                            gffexo.pop(max(gffexo.keys()))
                                    else:
                                        gffexo[max(gffexo.keys())]={min(gffexo[max(gffexo.keys())])+toshorten,
                                                                    max(gffexo[max(gffexo.keys())])}
                                    gffexo=sortgff(gffexo,sens)

# if we are working on another exon, no stop codon and no frame shift is allowed
                    elif "*" in exonaa or "#" in exonaa:
                        ispseudo2=True
                        break

# if ispseudo2 is true, then the hit is included inside pseudopossible
                if ispseudo2==True:
                    isalreadygood=False
                    if scaffold in exhits:
                        for h in exhits[scaffold]:
                            if exhits[scaffold][h]["debut"]<fin and exhits[scaffold][h]["fin"]>debut:
                                isalreadygood=True
                    if isalreadygood==False:
                        if scaffold not in pseudopossible:
                            pseudopossible[scaffold]={}
                        pseudopossible[scaffold][len(pseudopossible[scaffold])]={"debut":debut,
                                                "fin":fin,
                                                "score":int(rawscore)}
                    continue

# if ispseudo2 is false, the hit is added to exhits
# directly if the scafoold is not yet in exhits,
# otherwise it is compared with the other hits of the same scaffold
                if scaffold not in exhits:
                    exhits[scaffold]={}
                    exhits[scaffold][len(exhits[scaffold])]={"bestscore":rawscore,
                                                             "scaffold":scaffold,
                                                             "debut":int(debut),
                                                             "fin":int(fin),
                                                             "couverture":qlength/nlrs[query],
                                                             "gff":gffexo,
                                                             "rallonge":False,
                                                             "sens":sens}
                else:

# found is true if an overlapping hit is found in exhits, false otherwise
# the hit is then merged with the already-existing hit, and not added to exhits
                    found=False
                    
# this for loop is set to search for overlapping hit already in exhits,
# and if so to merge this existing hit with the currently analyzed hit
                    for h in exhits[scaffold].keys():
                        if exhits[scaffold][h]["debut"]<int(fin) and \
                            exhits[scaffold][h]["fin"]>int(debut) and \
                                exhits[scaffold][h]["sens"] == sens:
                            found=True
                            
# when merging two hits, we take the best value for "debut", "fin", "bestscore" and "couverture"
# "couverture" being the proportion of the query gene covered by the target gene
                            exhits[scaffold][h]["debut"]=min(exhits[scaffold][h]["debut"],int(debut))
                            exhits[scaffold][h]["fin"]=max(exhits[scaffold][h]["fin"],int(fin))
                            exhits[scaffold][h]["bestscore"]=max(exhits[scaffold][h]["bestscore"],rawscore)
                            exhits[scaffold][h]["couverture"]=max(exhits[scaffold][h]["couverture"],
                                                                  qlength/nlrs[query])
                            
# we take as gff the longest of the gffs of the two hits
                            if longueur(exhits[scaffold][h]["gff"])<longueur(gffexo):
                                exhits[scaffold][h]["gff"]=gffexo
                    
# if the current hit does not overlap with any already existing hit, we add it to exhits
                    if found==False:
                        exhits[scaffold][len(exhits[scaffold])]={"bestscore":rawscore,
                                                                 "scaffold":scaffold,
                                                                 "debut":int(debut),
                                                                 "fin":int(fin),
                                                                 "couverture":qlength/nlrs[query],
                                                                 "gff":gffexo,
                                                                 "rallonge":False,
                                                                 "sens":sens}

# we go there while we are reading the exonerate alignment between the query and the target
        elif alignt==True and totreat == True and ispseudo1 == False:
            
# if we reach "vulgar", it means the alignment reading is over
            if len(line.strip())>0 and line[0]!="#":
                elet=line.strip().split(":")
                if elet[0]=="vulgar":
                    alignt=False
                    continue

# cl count the line numbers inside the alignment
# cl%5=4 means we are reading an amino acid line of the target
# cl%5=0 means we are reading a nucleotides line of the target
# the sequence in the line is then named exonaa or exonnt
# it is split if a intro is included inside of it,
# and it is added to the corresponding list (cdsaa or cdsnt)
            cl+=1
            if cl%5==4:
                exonaa=line.strip().replace("{","").replace("}","").strip()
            if cl%5==0:
                exonnt=line.strip().split()[2].replace("{","").replace("}","").strip()
                pcxnt = exonnt.split(".")
                pcxaa = exonaa.split()
                for pc in pcxnt:
                    npc=pc
                    torem=[]
                    if npc != "":
                        for lettre in range(len(npc)):
                            if npc[lettre].islower() and pcxaa[0][lettre] in ["+","-"]:
                                torem.append(lettre)
                        for lettre in sorted(torem,reverse=True):
                            if lettre<len(npc)-1:
                                npc=npc[:lettre]+npc[lettre+1:]
                                pcxaa[0]=pcxaa[0][:lettre]+pcxaa[0][lettre+1:]
                            else:
                                npc=npc[:lettre]
                                pcxaa[0]=pcxaa[0][:lettre]
                        npc=npc.replace("-","")
                        pcxaa[0]=pcxaa[0].replace("-","").replace("+","")
                    if npc=="":
                        if cdsnt[-1]!="":
                            cdsnt.append("")
                            cdsaa.append("")
                        if len(pcxaa)>0:
                            if pcxaa[0]=="":
                                if len(pcxaa)>1:
                                    pcxaa=pcxaa[1:]
                                else:
                                    pcxaa=[]
                    else:
                        cdsnt[-1]=cdsnt[-1]+npc
                        ajouaa=pcxaa[0]
                        cdsaa[-1]=cdsaa[-1]+ajouaa
                        if len(pcxaa)>1:
                            pcxaa=pcxaa[1:]
                        else:
                            pcxaa=[]

# that condition is respected when when are reading the begining of a hit, before the alignment
# with general informations about the hit
        elif len(line.strip())>0 and line[0]!="#":
            elet=line.strip().split(":")
            if elet[0]=="Query":
                query="_".join(elet[1].strip().split("_")[:2])
            elif elet[0]=="Target":
                scaffold=elet[1].split()[0].strip()
            elif elet[0]=="Raw score":
                rawscore=elet[1].strip()
            elif elet[0]=="Query range":
                qlength=(int(elet[1].strip().split()[2])-int(elet[1].strip().split()[0])+1)*3
            elif elet[0]=="Target range":
                debut=min(int(elet[1].strip().split()[0]),int(elet[1].strip().split()[2]))
                fin=max(int(elet[1].strip().split()[0]),int(elet[1].strip().split()[2]))
                
# as "target range" is the last information before the alignment,
# here is tested if the hit overlaps with a NB,
# and if the NB is "bad, the hit is directly included in pseudonlr, and ispseudo1 is set to true
                totreat=False
                ispseudo1=False
                if scaffold in nb:
                    for i in nb[scaffold]:
                        if debut<nb[scaffold][i]["fin"] and fin>nb[scaffold][i]["debut"]:
                            totreat=True
                            if nb[scaffold][i]["status"]=="bad":
                                ispseudo1=True
                            indicenb=i
                if ispseudo1==True:
                    if scaffold not in pseudonlr:
                        pseudonlr[scaffold]={}
                        pseudonlr[scaffold][indicenb]={"debut":debut,
                                                "fin":fin,
                                                "score":int(rawscore),
                                                "genes":set()}
                    else:
                        alreadyin=False
                        if indicenb in pseudonlr[scaffold]:
                            alreadyin=True
                            pseudonlr[scaffold][indicenb]["debut"]=min(pseudonlr[scaffold][indicenb]["debut"],
                                                                       debut)
                            pseudonlr[scaffold][indicenb]["fin"]=min(pseudonlr[scaffold][indicenb]["fin"]
                                                                     ,fin)
                            pseudonlr[scaffold][indicenb]["score"]=max(pseudonlr[scaffold][indicenb]["score"],
                                                                int(rawscore))
                        else:
                            for i in pseudonlr[scaffold]:
                                if pseudonlr[scaffold][i]["debut"]<fin and pseudonlr[scaffold][i]["fin"]>debut:
                                    alreadyin=True
                                    pseudonlr[scaffold][i]["debut"]=min(pseudonlr[scaffold][i]["debut"],debut)
                                    pseudonlr[scaffold][i]["fin"]=min(pseudonlr[scaffold][i]["fin"],fin)
                                    pseudonlr[scaffold][i]["score"]=max(pseudonlr[scaffold][i]["score"],
                                                                        int(rawscore))
                                break
                        if alreadyin==False:
                            pseudonlr[scaffold][indicenb]={"debut":debut,
                                                    "fin":fin,
                                                    "score":int(rawscore),
                                                    "genes":set()}
                    
# that condition is respected when ispseudo1 is false
                elif totreat==True:
                    alignt=True
                    cl=0
                    cdsaa=[""]
                    cdsnt=[""]

# that condition is set to read the exonerate gff
# we read only the "cds" lines
            elif isgff==True:
                if line.strip().split()[2]=="cds":
                    sens=line.strip().split()[6]
                    gffexo[len(gffexo)]={int(line.strip().split()[3]),int(line.strip().split()[4])}
                                   
# once the exonerate output is fully read, we check if the hits in pseudopossible are overlapping with
# hits in exhits or in pseudonlr then the hit in pseudopossible is deleted,
# otherwise included in pseudonlr
for scaffold in pseudopossible:
    for i in pseudopossible[scaffold]:
        isnewpseudo=True
        if scaffold in exhits:
            for j in exhits[scaffold]:
                
# if the hit pseudopossible overlaps with a hit exhits, the hit is no more analyzed (and then deleted)
                if pseudopossible[scaffold][i]["debut"]<exhits[scaffold][j]["fin"] and \
                    pseudopossible[scaffold][i]["fin"]>exhits[scaffold][j]["debut"]:
                    isnewpseudo=False

# if it overlaps with a pseudonlr hit, it is merged with this pseudonlr hit
        if isnewpseudo==True:
            if scaffold in pseudonlr:
                for j in pseudonlr[scaffold]:
                    if pseudopossible[scaffold][i]["debut"]<pseudonlr[scaffold][j]["fin"] and \
                        pseudopossible[scaffold][i]["fin"]>pseudonlr[scaffold][j]["debut"]:
                        pseudonlr[scaffold][j]["debut"]=min(pseudopossible[scaffold][i]["debut"],
                                                                 pseudonlr[scaffold][j]["debut"])
                        pseudonlr[scaffold][j]["fin"]=min(pseudopossible[scaffold][i]["fin"],
                                                               pseudonlr[scaffold][j]["fin"])
                        pseudonlr[scaffold][j]["score"]=max(pseudopossible[scaffold][i]["score"],
                                                            pseudonlr[scaffold][j]["score"])
                        isnewpseudo=False

# if it does not overlap with any other hit, it is added to pseudonlr
        if isnewpseudo==True:
            if scaffold not in pseudonlr:
                pseudonlr[scaffold]={}
                pseudonlr[scaffold][0]={
                    "debut":pseudopossible[scaffold][i]["debut"],
                    "fin":pseudopossible[scaffold][i]["fin"],
                    "score":pseudopossible[scaffold][i]["score"],
                    "genes":set()}
            else:
                pseudonlr[scaffold][max(pseudonlr[scaffold].keys())+1]={
                    "debut":pseudopossible[scaffold][i]["debut"],
                    "fin":pseudopossible[scaffold][i]["fin"],
                    "score":pseudopossible[scaffold][i]["score"],
                    "genes":set()}

del pseudopossible

with open("./readable_"+ind+"_correc_exo.txt",'a') as readable:    
    readable.write("#whole-NLRs exonerate output analyzed\n")
    readable.write("#time is "+str(timedelta(seconds=(timer()-t0)))+"\n#\n")
    readable.write("#\n#merging of the pseudo-NLR genes\n")

###############################################################################
######## 4. merging the pseudo-NLR genes
###############################################################################

# here the pseudo-NLR genes are merged if they overlap
# and the bast informations are kept when merging the two pseudo-NLRs
for scaffold in pseudonlr:
    toremove=set()
    for i in pseudonlr[scaffold]:
        if i in toremove:continue
        isok=False
        while isok!=True:
            isok=True
            for j in pseudonlr[scaffold]:
                if j==i or j in toremove:continue
                if pseudonlr[scaffold][i]["debut"]<pseudonlr[scaffold][j]["fin"] and \
                    pseudonlr[scaffold][i]["fin"]>pseudonlr[scaffold][j]["debut"]:
                    pseudonlr[scaffold][i]["debut"]=min(pseudonlr[scaffold][i]["debut"],
                                                             pseudonlr[scaffold][j]["debut"])
                    pseudonlr[scaffold][i]["fin"]=min(pseudonlr[scaffold][i]["fin"],
                                                           pseudonlr[scaffold][j]["fin"])
                    pseudonlr[scaffold][i]["score"]=max(pseudonlr[scaffold][i]["score"],
                                                        pseudonlr[scaffold][j]["score"])
                    toremove.add(j)
                    isok=False
    for i in toremove:
        del pseudonlr[scaffold][i]


with open("./readable_"+ind+"_correc_exo.txt",'a') as readable:    
    readable.write("#pseudo-NLR genes have been merged\n")
    readable.write("#time is "+str(timedelta(seconds=(timer()-t0)))+"\n#\n")
    readable.write("#\n#analyse of the gffs of existing genes\n")

###############################################################################
######## 5. looking for overlapping genes and their gff
###############################################################################

# here we search for overlapping genes with the hits

# genes is the dico per scaffold of the genes overlapping with hits in exhits 
genes={}

# we read a first time the pre-existing gff of the individual,
# this time we are only interested by the genes and their limits
if args.gff_file:
    with open(args.gff_file,"r") as gff_file:
        for gline in gff_file:
            if gline.strip().split()[2]=="gene":
                scaffold=gline.strip().split()[0]
                
    # we are first interested by genes overlapping with a hit in exhits
                if scaffold in exhits.keys():
                    for i in exhits[scaffold].keys():
                        if int(gline.strip().split()[3])<exhits[scaffold][i]["fin"] and \
                            int(gline.strip().split()[4])>exhits[scaffold][i]["debut"] and \
                                gline.strip().split()[6] == exhits[scaffold][i]["sens"]: 
                            if scaffold not in genes.keys():
                                genes[scaffold]={}
                            if gline.strip().split()[8].replace("ID=","") not in genes[scaffold].keys():
                                genes[scaffold][gline.strip().split()[8].replace("ID=","")]={
                                    "debut":int(gline.strip().split()[3]),
                                    "fin":int(gline.strip().split()[4]),
                                    "status":"good",
                                    "sens":gline.strip().split()[6]}
                                if ind+"_"+gline.strip().split()[8].replace("ID=","") in nlr_new:
                                    genes[scaffold][gline.strip().split()[8].replace("ID=","")]["status"]="bad"
                        elif ind+"_"+gline.strip().split()[8].replace("ID=","") in nlrlist_ind: 
                            if scaffold not in genes.keys():
                                genes[scaffold]={}
                            if gline.strip().split()[8].replace("ID=","") not in genes[scaffold].keys():
                                genes[scaffold][gline.strip().split()[8].replace("ID=","")]={
                                    "debut":int(gline.strip().split()[3]),
                                    "fin":int(gline.strip().split()[4]),
                                    "status":"good",
                                    "sens":gline.strip().split()[6]}
                                if ind+"_"+gline.strip().split()[8].replace("ID=","") in nlr_new:
                                    genes[scaffold][gline.strip().split()[8].replace("ID=","")]["status"]="bad"
    
    # we are secondarily interested by genes overlapping with a hit in pseudonlr
                if scaffold in pseudonlr.keys():
                    for i in pseudonlr[scaffold].keys():
                        if int(gline.strip().split()[3])<pseudonlr[scaffold][i]["fin"] and \
                            int(gline.strip().split()[4])>pseudonlr[scaffold][i]["debut"]: 
                            if scaffold not in genes.keys():
                                genes[scaffold]={}
                            if gline.strip().split()[8].replace("ID=","") not in genes[scaffold].keys():
                                genes[scaffold][gline.strip().split()[8].replace("ID=","")]={
                                    "debut":int(gline.strip().split()[3]),
                                    "fin":int(gline.strip().split()[4]),
                                    "status":"good",
                                    "sens":gline.strip().split()[6]}
                                if ind+"_"+gline.strip().split()[8].replace("ID=","") in nlr_new:
                                    genes[scaffold][gline.strip().split()[8].replace("ID=","")]["status"]="bad"

# oldgff is the dico of the gffs 
oldgff={}

# we read a second time the pre-existing gff of the individual,
# this time we extract the gff of the genes in the dico genes
if args.gff_file:
    with open(args.gff_file,"r") as gff_file:
        for gline in gff_file:
    
    # we take only the CDS in consideration
            if gline.strip().split()[2]=="CDS":
                
    # parent is the name of the gene "parent" of the CDS, possibly in the genes dico
                parent={s.split("_")[0].replace("Parent=","") for s in gline.strip().split()[8].split(";") \
                        if "Parent=" in s}
                parent=parent.pop()
                scaffold=gline.strip().split()[0]
                if scaffold in genes.keys():
                    if parent in genes[scaffold].keys():
                        
    # we add the CDS to oldgff if the parent is in the genes dico
                        if parent not in oldgff:
                            oldgff[parent]={}
                        oldgff[parent][len(oldgff[parent])]={int(gline.strip().split()[3]),
                                                             int(gline.strip().split()[4])}

with open("./readable_"+ind+"_correc_exo.txt",'a') as readable:    
    readable.write("#existing genes searched\n")
    readable.write("#time is "+str(timedelta(seconds=(timer()-t0)))+"\n#\n")
    readable.write("#\n#checking the status of those existing genes\n")

###############################################################################
######## 6. checking overlapping genes status : good or bad
###############################################################################

# for a gene to be considered good, it needs to overlap with a NB
# and it needs to not contain any non-terminal stop codon
for scaffold in genes:
    for record in SeqIO.parse(args.assembly_file,"fasta"):
        if str(record.id)==scaffold:
            for g in genes[scaffold]:

# we first check for the presence of a NB, good or bad, overlapping with the given gene
                nbpresent=False
                if scaffold in nb:
                    for i in nb[scaffold]:
                        if genes[scaffold][g]["debut"]<nb[scaffold][i]["fin"] and \
                            genes[scaffold][g]["fin"]>nb[scaffold][i]["debut"]:
                            nbpresent=True
                if nbpresent==False:
                    genes[scaffold][g]["status"]="bad"
                    continue

# in a second time, we check for the integrity of the sequence
                wholeseq=""
                gf=sortgff(oldgff[g],genes[scaffold][g]["sens"])
                for i in sorted(list(gf.keys())):
                    if len(gf[i])>1:
                        sq=record.seq[min(gf[i])-1:max(gf[i])]
                    else:
                        sq=Seq(record.seq[min(gf[i])-1])
                    if genes[scaffold][g]["sens"]=="-":
                        sq=sq.reverse_complement()
                    wholeseq=wholeseq+str(sq)
                prot=str(Seq(wholeseq).translate())
                if "*" in prot[:-1]:
                    genes[scaffold][g]["status"]="bad"
                if len(wholeseq)%3!=0:
                    if genes[scaffold][g]["sens"]=="-" and int(min(gf[max(gf.keys())]))!=1:
                        genes[scaffold][g]["status"]="bad"
                    elif genes[scaffold][g]["sens"]=="+" and int(max(gf[max(gf.keys())])) != len(record.seq):
                        genes[scaffold][g]["status"]="bad"

with open("./readable_"+ind+"_correc_exo.txt",'a') as readable:    
    readable.write("#status of pre-existing genes identified\n")
    readable.write("#time is "+str(timedelta(seconds=(timer()-t0)))+"\n#\n")
    readable.write("length of exhits :"+str(sum([len(exhits[sca]) for sca in exhits]))+"\n")
    readable.write("length of genes :"+str(len(oldgff))+"\n")
    readable.write("#\n#merging of genes and exonerate hits\n")

###############################################################################
######## 7. making intersection between genes and good hits
###############################################################################              

# hits contains all the genes overlapping with a hit, and a consensus of both of them
# hits not overlapping with a gene are storred per scaffold in the "unknown" sub-dico
hits={"unknown":{}}

# scaflist is the set of all the scaffold in exhits
scaflist=set(exhits.keys()).copy()

for scaffold in scaflist:
    
# keylist is the set of all the keys corresponding to hits in exhits, for a given scaffold
    keylist=set(exhits[scaffold].keys()).copy()
    for h in keylist:
        
# found is false if at the end of the turn the hit h is to be added at "unknown",
# true if it associated with a pre-existing gene
        found=False
        if scaffold in genes:
            for g in genes[scaffold]:
                
# we check if the considered gene overlaps with the considered hit
                if genes[scaffold][g]["debut"]<exhits[scaffold][h]["fin"] and \
                    genes[scaffold][g]["fin"]>exhits[scaffold][h]["debut"] and \
                        exhits[scaffold][h]["sens"] == genes[scaffold][g]["sens"]:                    
                    found=True
                    
# there are two posssible cases : the gene is already in hits, or not
# if not, it is added to hits
                    if ind+"_"+g not in hits:
                        
# the gene is first added with the informations of the exonerate hit
# if the gene status is "good", hit and gene are merged, priorizing the hit over the gene
                        hits[ind+"_"+g]=exhits[scaffold][h]
                        if genes[scaffold][g]["status"]=="good":
                            if genes[scaffold][g]["debut"]>exhits[scaffold][h]["debut"] and \
                                genes[scaffold][g]["fin"]<exhits[scaffold][h]["fin"]:
                                hits[ind+"_"+g]["rallonge"]=True
                            hits[ind+"_"+g]["debut"]=min(exhits[scaffold][h]["debut"],
                                                         genes[scaffold][g]["debut"])
                            hits[ind+"_"+g]["fin"]=max(exhits[scaffold][h]["fin"],
                                                       genes[scaffold][g]["fin"])
                            
# there we call the mergegff function to merge the gff of the gene (oldgff) and the exonerate gff
                            for record in SeqIO.parse(args.assembly_file,"fasta"):
                                if str(record.id)==scaffold:
                                    hits[ind+"_"+g]["gff"]=mergegff(exhits[scaffold][h]["gff"],
                                                                       oldgff[g],
                                                                       exhits[scaffold][h]["sens"],
                                                                       record)

# if the gene is already in hits, we simply merge the considered exonerate hit with the hit in hits
# we then priorize the longest of both, calling the longueur function for that
                    else:
                        if hits[ind+"_"+g]["bestscore"]<exhits[scaffold][h]["bestscore"]:
                            hits[ind+"_"+g]["bestscore"]=exhits[scaffold][h]["bestscore"]
                            hits[ind+"_"+g]["couverture"]=exhits[scaffold][h]["couverture"]
                        if hits[ind+"_"+g]["debut"]>exhits[scaffold][h]["debut"]:
                            hits[ind+"_"+g]["debut"]=exhits[scaffold][h]["debut"]
                            hits[ind+"_"+g]["rallonge"]=True
                        if hits[ind+"_"+g]["fin"]<exhits[scaffold][h]["fin"]:
                            hits[ind+"_"+g]["fin"]=exhits[scaffold][h]["fin"]
                            hits[ind+"_"+g]["rallonge"]=True
                        if longueur(hits[ind+"_"+g]["gff"])>longueur(exhits[scaffold][h]["gff"]):
                            for record in SeqIO.parse(args.assembly_file,"fasta"):
                                if str(record.id)==scaffold:
                                    hits[ind+"_"+g]["gff"]=mergegff(hits[ind+"_"+g]["gff"],
                                                                    exhits[scaffold][h]["gff"],
                                                                    hits[ind+"_"+g]["sens"],
                                                                    record)
                        else:
                            for record in SeqIO.parse(args.assembly_file,"fasta"):
                                if str(record.id)==scaffold:
                                    hits[ind+"_"+g]["gff"]=mergegff(exhits[scaffold][h]["gff"],
                                                                    hits[ind+"_"+g]["gff"],
                                                                    hits[ind+"_"+g]["sens"],
                                                                    record)
                            
# if the hit is not found as overlapping with an already existing gene
# we add it to "unknown"
        if found==False:
            if scaffold not in hits["unknown"]:
                hits["unknown"][scaffold]={}
            hits["unknown"][scaffold][len(hits["unknown"][scaffold])]=exhits[scaffold][h]
        
# we remove the hit from exhits to gain some memory
        del exhits[scaffold][h]

with open("./readable_"+ind+"_correc_exo.txt",'a') as readable:    
    readable.write("#genes and good exonerate hits merged\n")
    readable.write("#time is "+str(timedelta(seconds=(timer()-t0)))+"\n#\n")
    readable.write("#\n#merging good 'unknown' hits now\n")

###############################################################################
######## 8. merging of good unknwon hits between then
###############################################################################            

for scaffold in hits["unknown"]:
    
# first thing is that we sort all the unknown hits of a given scaffold, based on their value of "debut"
    newunknown={}
    for ing,iu in zip(range(len(hits["unknown"][scaffold])),
                   [iun[0] for iun in sorted(hits["unknown"][scaffold].items(),
                                             key=lambda item: item[1]["debut"])]):
        newunknown[ing]=hits["unknown"][scaffold][iu]
    hits["unknown"][scaffold]=newunknown
    
# positions is the set of all the keys of unknwon hits of the given scaffold
    positions=set(hits["unknown"][scaffold].keys())

# newlist is the dico of all the unknown hits on the scaffold that will be keep finally
    newlist={}
    
    for i in range(len(hits["unknown"][scaffold])):
        
# this condition is set for working only on hits we haven't work yet on
        if i in positions:
            positions.remove(i)
            
            # we initiate values based on the considered unknwon hit
            debut=hits["unknown"][scaffold][i]["debut"]
            fin=hits["unknown"][scaffold][i]["fin"]
            score=hits["unknown"][scaffold][i]["bestscore"]
            couverture=hits["unknown"][scaffold][i]["couverture"]
            gff=hits["unknown"][scaffold][i]["gff"]
            sens=hits["unknown"][scaffold][i]["sens"]
            
# we compare each hit to the next ones on the scaffold
# we keep the longest values for start, end and score
# we keep only the longest gff (we do not merge them)
            for k in range(i+1,len(hits["unknown"][scaffold])):
                if k in positions:
                    if hits["unknown"][scaffold][k]["debut"]<fin and hits[
                            "unknown"][scaffold][k]["fin"]>debut and hits[
                                "unknown"][scaffold][k]["sens"]==sens:
                        positions.remove(k)
                        if hits["unknown"][scaffold][k]["debut"]<debut:
                            debut=hits["unknown"][scaffold][k]["debut"]
                        if hits["unknown"][scaffold][k]["fin"]>fin:
                            fin=hits["unknown"][scaffold][k]["fin"]
                        if hits["unknown"][scaffold][k]["bestscore"]>score:
                            score=hits["unknown"][scaffold][k]["bestscore"]
                            couverture=hits["unknown"][scaffold][k]["couverture"]
                        if longueur(hits["unknown"][scaffold][k]["gff"])>longueur(gff):
                            gff=hits["unknown"][scaffold][k]["gff"]
            
# we add the new "reconstructed" hit to newlist
            newlist[len(newlist)]={"bestscore":score,
                                    "debut":debut,
                                    "fin":fin,
                                    "gff":gff,
                                    "couverture":couverture,
                                    "sens":sens}

# at the end, newlist take place as hits["unknown"][scaffold]
    hits["unknown"][scaffold]=newlist

with open("./readable_"+ind+"_correc_exo.txt",'a') as readable:    
    readable.write("#'unknown' hits merged\n")
    readable.write("#time is "+str(timedelta(seconds=(timer()-t0)))+"\n#\n")
    readable.write("#comparison of 'unknown' hits and gene hits now\n")
  
###############################################################################
######## 9. merging of good unknwon hits with good gene hits
###############################################################################            

# "hitslist" is a set of all the genes in hits
hitslist=set(hits.keys()).copy()

for h in hitslist:
    if h=="unknown":continue
    scaffold=hits[h]["scaffold"]
    
# we compare each gene with all the unknown hits of its scaffold
# if a gene can be prolongated by an unknown hit, its value "rallonge" is set to true
# and the gffs are merged, priorizing the gene over the unknwon hit
    if scaffold in hits["unknown"]:

# toremove is a set containing all the unknown hits to remove because overlapping with a gene hit
        toremove=set()
        
        for i in hits["unknown"][scaffold].keys():
            if hits["unknown"][scaffold][i]["debut"]<hits[h]["fin"] and \
                hits["unknown"][scaffold][i]["fin"]>hits[h]["debut"] and \
                    hits["unknown"][scaffold][i]["sens"] == hits[h]["sens"]:
                if hits["unknown"][scaffold][i]["fin"]>hits[h]["fin"]:
                    hits[h]["fin"]=hits["unknown"][scaffold][i]["fin"]
                    hits[h]["rallonge"]=True
                if hits["unknown"][scaffold][i]["debut"]<hits[h]["debut"]:
                    hits[h]["debut"]=hits["unknown"][scaffold][i]["debut"]
                    hits[h]["rallonge"]=True
                for record in SeqIO.parse(args.assembly_file,"fasta"):
                    if str(record.id)==scaffold:
                        hits[h]["gff"]=mergegff(hits[h]["gff"],
                                                hits["unknown"][scaffold][i]["gff"],
                                                hits[h]["sens"],
                                                record)
                toremove.add(i)
        for i in toremove:
            hits["unknown"][scaffold].pop(i)

with open("./readable_"+ind+"_correc_exo.txt",'a') as readable:                    
    readable.write("#hits 'gene' and 'unknown' merged\n")
    readable.write("#time is "+str(timedelta(seconds=(timer()-t0)))+"\n#\n")
    readable.write("#comparison of gene hits between them now\n")

###############################################################################
######## 10. merging of good gene hits between them
###############################################################################            
    
toremove=set()

for h in hitslist:
    if h=="unknown":continue

# this condition is set because we work on each gene only once
    if h in toremove:continue
    
# we initialize the values with those of the studied gene
    start=hits[h]["debut"]
    end=hits[h]["fin"]
    scaffold=hits[h]["scaffold"]
    score=hits[h]["bestscore"]
    couverture=hits[h]["couverture"]
    gff=hits[h]["gff"]
    rallonge=hits[h]["rallonge"]
    sens=hits[h]["sens"]

# newadd is false if after a turn of the while loop we found no overlapping gene, true otherwise
    newadd=True

# toadd correspond to the studied gene "h" and all the genes overlapping with it
    toadd=set()
    while newadd==True:
        newadd=False
        for k in hitslist:
            if k=="unknown":continue

# the gene hit k is compared to the gene hit h, to do so, it must be :
# different of h, not to add to h, not to remove (already worked on), on the same scaffold as h
# and of course it must overlap with h and be on the same strand
            if k != h and k not in toadd.union(toremove) and scaffold==hits[k]["scaffold"]:
                if hits[k]["debut"]<end and hits[k]["fin"]>start and hits[k]["sens"] == sens:
                    toadd.add(k)
                    if h not in toadd:
                        toadd.add(h)
                    newadd=True
                    start=min(start,hits[k]["debut"])
                    end=max(end,hits[k]["fin"])
                    score=max(score,hits[k]["bestscore"])
                    couverture=max(couverture,hits[k]["couverture"])
                    for record in SeqIO.parse(args.assembly_file,"fasta"):
                        if str(record.id)==scaffold:
                            if longueur(gff)>=longueur(hits[k]["gff"]):
                                gff=mergegff(gff,hits[k]["gff"],sens,record)
                            else:
                                gff=mergegff(hits[k]["gff"],gff,sens,record)
                    rallonge=True

# this condition is satisfied only if any other gene hit has been found as overlapping with the hit gene h
# a new hit, intersection of all the overlapping hit genes is then created,
# and the single gene hits are deleted
    if len(toadd)>0:
        toremove=toremove.union(toadd)
        hits["+".join(toadd)]={"bestscore":score,
                                "scaffold":scaffold,
                                "debut":start,
                                "fin":end,
                                "couverture":couverture,
                                "gff":gff,
                                "rallonge":True,
                                "sens":sens}
for i in toremove:
    hits.pop(i)

with open("./readable_"+ind+"_correc_exo.txt",'a') as readable:
    readable.write("#hits 'gene' merged\n")
    readable.write("#time is "+str(timedelta(seconds=(timer()-t0)))+"\n#\n")
    readable.write("#comparison of gene hits between them now\n")

###############################################################################
######## 11. looking for bad hits being known genes
###############################################################################              

# we look for genes in genes overlapping with bad hits in pseudonlr
for scaffold in pseudonlr:
    if scaffold in genes:
        for i in pseudonlr[scaffold]:
            for g in genes[scaffold]:

# we considered that a pseudo-NLR overlapping with a gene being a good NLR means that
# this gene was a chimera containing both good and bad NLRs
# the name of the gene is then only kept for the good NLR
# (but this situation has very low probability to happen)
# (it must then not be problematic since the gene should be considered as bad given the situation)
                if ind+"_"+g in hits: continue
                if genes[scaffold][g]["debut"]<pseudonlr[scaffold][i]["fin"] and \
                    genes[scaffold][g]["fin"]>pseudonlr[scaffold][i]["debut"]:
                    pseudonlr[scaffold][i]["genes"].add(ind+"_"+g)
                    pseudonlr[scaffold][i]["debut"]=min(pseudonlr[scaffold][i]["debut"],
                                                             genes[scaffold][g]["debut"])
                    pseudonlr[scaffold][i]["fin"]=min(pseudonlr[scaffold][i]["fin"],
                                                           genes[scaffold][g]["fin"])

# to finish with pseudonlr, we merge the rare pseu-NLR genes that could still be overlapping
# (just to be sure)
for scaffold in pseudonlr:
    toremove=set()
    for i in pseudonlr[scaffold]:
        if i in toremove:continue
        isok=False
        while isok!=True:
            isok=True
            for j in pseudonlr[scaffold]:
                if j==i or j in toremove:continue
                if pseudonlr[scaffold][i]["debut"]<pseudonlr[scaffold][j]["fin"] and \
                    pseudonlr[scaffold][i]["fin"]>pseudonlr[scaffold][j]["debut"]:
                    pseudonlr[scaffold][i]["debut"]=min(pseudonlr[scaffold][i]["debut"],
                                                             pseudonlr[scaffold][j]["debut"])
                    pseudonlr[scaffold][i]["fin"]=min(pseudonlr[scaffold][i]["fin"],
                                                           pseudonlr[scaffold][j]["fin"])
                    pseudonlr[scaffold][i]["score"]=max(pseudonlr[scaffold][i]["score"],
                                                        pseudonlr[scaffold][j]["score"])
                    pseudonlr[scaffold][i]["genes"]=pseudonlr[scaffold][i]["genes"].union(
                        pseudonlr[scaffold][j]["genes"])
                    toremove.add(j)
                    isok=False
    for i in toremove:
        del pseudonlr[scaffold][i]

nb_pseudo=0
for scaffold in pseudonlr:
    nb_pseudo+=len(pseudonlr[scaffold])



with open("./readable_"+ind+"_correc_exo.txt",'a') as readable:    
    readable.write("#genes and bad exonerate hits merged\n")
    readable.write("#time is "+str(timedelta(seconds=(timer()-t0)))+"\n#\n")
    readable.write("#\n#checking previously-known NLRs\n")

###############################################################################
######## 12. checking previously-known NLRs
###############################################################################      

# nlrnonb is a set of all the previous-nlr not identified with exonerate (because no NB identified)
nlrnonb=set()
newgenelist=set(hits.keys()).copy()
newgenelist.remove("unknown")

for oldnlr in nlrlist_ind:
    nbpresent=False
    israllonge=True
    ispseudo=False
    hitslist=set(hits.keys()).copy()
    
# we first check if the old NLR is in hits
    for k in hitslist:
        if k=="unknown":continue
        if "+" in k:
            hl=set(k.split("+"))
        else:
            hl={k}
        for h in hl:
            if h == oldnlr:
                nlrlist_prolong.add(k)
                if k in newgenelist:
                    newgenelist.remove(k)
                nbpresent=True
                scaffold=hits[k]["scaffold"]
                if hits[k]["rallonge"]==False and genes[scaffold][oldnlr.split("_")[1]]["status"]=="good":
                    israllonge=False
                    if k in nlrlist_prolong:
                        nlrlist_prolong.remove(k)
                    hits.pop(k)
                    break

# then we check if the old NLR is in pseudonlr
    if nbpresent==False:
        for sc in pseudonlr:
            pseudolist=set(pseudonlr[sc]).copy()
            for k in pseudolist:
                if oldnlr in pseudonlr[sc][k]["genes"]:
                    nbpresent==True
                    ispseudo=True
                    if genes[sc][oldnlr.split("_")[1]]["status"]=="good":
                        pseudonlr[sc].pop(k)
                        
# nbpresent is false if the gene is not in hits nor in pseudonlr
    if nbpresent==False:
        nlrnonb.add(oldnlr)
        
# israllonge is false if the gene is in hits and not prolongated by an exonerate hit
    elif israllonge==False:
        nlrlist_notouch[oldnlr]=scaffold

# israllonge is false if the gene is in hits and not prolongated by an exonerate hit
    elif ispseudo==True:
        nlrlist_pseudo.add(oldnlr)

# we do the same with sacer genes
for sace in sacers.keys():
    isfound=False
    hitslist=set(hits.keys()).copy()
    for k in hitslist:
        if k=="unknown":continue
        if "+" in k:
            hl=set(k.split("+"))
        else:
            hl={k}
        for h in hl:
            if h == sace:
                isfound=True
                if k not in nlrlist_prolong:
                    sacerlist_kept.add(k)
                if k in newgenelist:
                    newgenelist.remove(k)
                break
    if isfound==False:
        sacer_notkept.add(sace)

with open("./readable_"+ind+"_correc_exo.txt",'a') as readable: 
    readable.write("#previously-known NLRs checked\n")
    readable.write("#time is "+str(timedelta(seconds=(timer()-t0)))+"\n#\n")
    readable.write("#\n#extraction of the sequences, and check of theit length\n")
   
###############################################################################
######## 13. extraction of the sequences
###############################################################################

# first we look for the sequences of the nlrs which have not been affected by the correction
# those two dicos are the fastas of these genes, repsectively in nucleotides and in amino acids
fastant_nlrnotouch={}
fastaaa_nlrnotouch={}
nlrlist=set(nlrlist_notouch.keys()).copy()
for g in nlrlist:

# nlrsq is the sequence of the given gene in nt
    nlrsq=""
    gf=sortgff(oldgff[g.split("_")[1]],genes[nlrlist_notouch[g]][g.split("_")[1]]["sens"])

# this part is used to build the sequence nlrsq from the gff gf
    for i in range(len(gf)):
        for record in SeqIO.parse(args.assembly_file,"fasta"):
            if str(record.id)==nlrlist_notouch[g]:
                if len(gf[i])>1:
                    sq=record.seq[min(gf[i])-1:max(gf[i])]
                else:
                    sq=Seq(record.seq[min(gf[i])-1])
                if genes[nlrlist_notouch[g]][g.split("_")[1]]["sens"]=="-":
                    sq=sq.reverse_complement()
                nlrsq=nlrsq+str(sq)

# prot is the amino acid sequence (translated from nlrsq)
    prot=str(Seq(nlrsq).translate())

# prot is not supposed to contain any stop codon, but at its extremity
# if so, there was a bug previously
    if "*" in prot[:-1]:
        print("problem : there is a stop codon in the sequence of gene not touch "+g)
        print("prot sequence : "+prot)
        print("nt sequence : "+nlrsq)
        print("gff : ")
        print(gf)
        os.kill()

# nlrsq is supposed to have a length multiple of 3
# if not, there as bug previously, or an old gff is wrong
    if len(nlrsq)%3!=0:
        if genes[nlrlist_notouch[g]][g.split("_")[1]]["sens"]=="-" and int(min(gf[max(gf.keys())]))!=1:
            print("problem : the sequence of gene not touch "+g+" is not a multiple of 3")
            print("prot sequence : "+prot)
            print("nt sequence : "+nlrsq)
            print("gff : ")
            print(gf)
            os.kill()
        elif genes[nlrlist_notouch[g]][g.split("_")[1]]["sens"]=="+":
            for record in SeqIO.parse(args.assembly_file,"fasta"):
                if str(record.id)==nlrlist_notouch[g]:
                    if int(max(gf[max(gf.keys())])) != len(record.seq):
                        print("problem : the sequence of gene not touch "+g+" is not a multiple of 3")
                        print("prot sequence : "+prot)
                        print("nt sequence : "+nlrsq)
                        print("gff : ")
                        print(gf)
                        os.kill()

# if the sequence does not end with a stop codon, in this while loop we had successively
# the following codon to the end of athe sequence, until this codon is a stop codon
# or until we reach the end of the scaffold
# we modify the gff consequently
    if prot[-1]!="*":
        lastex=max(gf.keys())
        newcodon=""
        reachtheend=False
        while str(Seq(newcodon).translate())!="*" and reachtheend==False:
            for record in SeqIO.parse(args.assembly_file,"fasta"):
                if str(record.id)==nlrlist_notouch[g]:
                    if genes[nlrlist_notouch[g]][g.split("_")[1]]["sens"]=="+":
                        if len(record.seq)<max(gf[lastex])+3:
                            reachtheend=True
                            break
                        newcodon=str(record.seq[max(gf[lastex]):max(gf[lastex])+3])
                        gf[lastex]={min(gf[lastex]),max(gf[lastex])+3}
                        nlrsq=nlrsq+newcodon
                        prot=prot+str(Seq(newcodon).translate())
                    elif genes[nlrlist_notouch[g]][g.split("_")[1]]["sens"]=="-":
                        if min(gf[lastex])-4<0:
                            reachtheend=True
                            break
                        newcodon=str(record.seq[min(gf[lastex])-4:min(gf[lastex])-1].reverse_complement())
                        gf[lastex]={min(gf[lastex])-3,max(gf[lastex])}
                        nlrsq=nlrsq+newcodon
                        prot=prot+str(Seq(newcodon).translate())


# this condition is used to be sure that the protein begins with a methionin amino acid
    if prot[0]!="M":
        endloop=False
        while endloop==False:
            for record in SeqIO.parse(args.assembly_file,"fasta"):
                if str(record.id)==nlrlist_notouch[g]:
                    if genes[nlrlist_notouch[g]][g.split("_")[1]]["sens"]=="+":
                        if min(gf[0])-3<1:
                            endloop=True
                            break
                        newcodon=str(record.seq[min(gf[0])-4:min(gf[0])-1])
                        if str(Seq(newcodon).translate())=="*":
                            endloop=True
                            break
                        else:
                            gf[0]={min(gf[0])-3,max(gf[0])}
                            nlrsq=newcodon+nlrsq
                            prot=str(Seq(newcodon).translate())+prot
                            if str(Seq(newcodon).translate())=="M":
                                endloop=True
                                break
                    elif genes[nlrlist_notouch[g]][g.split("_")[1]]["sens"]=="-":
                        if max(gf[0])+3>len(record.seq):
                            endloop=True
                            break
                        newcodon=str(record.seq[max(gf[0]):max(gf[0])+3])
                        if str(Seq(newcodon).translate().reverse_complement())=="*":
                            endloop=True
                            break
                        else:
                            gf[0]={min(gf[0]),max(gf[0])+3}
                            nlrsq=newcodon+nlrsq
                            prot=str(Seq(newcodon).translate().reverse_complement())+prot
                            if str(Seq(newcodon).translate().reverse_complement())=="M":
                                endloop=True
                                break
        if prot[0]!="M":
            if prot.find("M")==-1:
                nlrlist_notouch.pop(g)
                continue
            else:
                toremove=prot.find("M")*3
                nlrsq=nlrsq[toremove:]
                prot=prot[int(toremove/3):]
                if genes[nlrlist_notouch[g]][g.split("_")[1]]["sens"]=="+":
                    indices=sorted(list(gf.keys())).copy()
                    for i in indices:
                        if toremove>=(max(gf[i])-min(gf[i])+1):
                            toremove=toremove-(max(gf[i])-min(gf[i])+1)
                            gf.pop(i)
                            if toremove==0:
                                gf=sortgff(gf,genes[nlrlist_notouch[g]][g.split("_")[1]]["sens"])
                                break
                        else:
                            gf[i]={min(gf[i])+toremove,max(gf[i])}
                            break
                elif genes[nlrlist_notouch[g]][g.split("_")[1]]["sens"]=="-":
                    indices=sorted(list(gf.keys())).copy()
                    for i in indices:
                        if toremove>=(max(gf[i])-min(gf[i])+1):
                            toremove=toremove-(max(gf[i])-min(gf[i])+1)
                            gf.pop(i)
                            if toremove==0:
                                gf=sortgff(gf,genes[nlrlist_notouch[g]][g.split("_")[1]]["sens"])
                                break
                        else:
                            gf[i]={min(gf[i]),max(gf[i])-toremove}
                            break

# this final condition is used to filtrate the genes too short to be complete NLRs
    l=len(nlrsq)
    if l<int(args.minsize_nlr):
        nlrlist_notouch.pop(g)
    else:
        oldgff[g]=gf
        fastant_nlrnotouch[g]=nlrsq
        fastaaa_nlrnotouch[g]=prot

# we do the same we did just before, but for genes in hits
fastant_hits={}
fastaaa_hits={}
nlrlist=nlrlist_prolong.union(sacerlist_kept).union(newgenelist).copy()
for g in nlrlist:
    nlrsq=""
    gf=sortgff(hits[g]["gff"],hits[g]["sens"])
    for i in range(len(gf)):
        for record in SeqIO.parse(args.assembly_file,"fasta"):
            if str(record.id)==hits[g]["scaffold"]:
                if len(gf[i])>1:
                    sq=record.seq[min(gf[i])-1:max(gf[i])]
                else:
                    sq=Seq(record.seq[min(gf[i])-1])
                if hits[g]["sens"]=="-":
                    sq=sq.reverse_complement()
                nlrsq=nlrsq+str(sq)
    prot=str(Seq(nlrsq).translate())
    if "*" in prot[:-1]:
        print("problem : there is a stop codon in the sequence of hit "+g)
        print("prot sequence : "+prot)
        print("nt sequence : "+nlrsq)
        print("gff : ")
        print(gf)
        os.kill()

# a difference with before is that here, if a gene length is a not a multiple of 3,
# here we just add the one or two missing nucleotides
    if len(nlrsq)%3!=0:
        lastex=max(gf.keys())
        missing=3-(len(nlrsq)%3)
        for record in SeqIO.parse(args.assembly_file,"fasta"):
            if str(record.id)==hits[g]["scaffold"]:
                if hits[g]["sens"]=="+":
                    if len(record.seq)<max(gf[lastex])+3:
                        break
                    newcodon=str(record.seq[max(gf[lastex]):max(gf[lastex])+missing])
                    gf[lastex]={min(gf[lastex]),max(gf[lastex])+missing}
                    nlrsq=nlrsq+newcodon
                    prot=prot+str(Seq(newcodon).translate())
                elif hits[g]["sens"]=="-":
                    if min(gf[lastex])-4<0:
                        break
                    newcodon=str(record.seq[min(gf[lastex])-1-missing:min(gf[lastex])-1].reverse_complement())
                    gf[lastex]={min(gf[lastex])-missing,max(gf[lastex])}
                    nlrsq=nlrsq+newcodon
                    prot=prot+str(Seq(newcodon).translate())
    if prot[-1]!="*":
        lastex=max(gf.keys())
        newcodon=""
        reachtheend=False
        while str(Seq(newcodon).translate())!="*" and reachtheend==False:
            for record in SeqIO.parse(args.assembly_file,"fasta"):
                if str(record.id)==hits[g]["scaffold"]:
                    if hits[g]["sens"]=="+":
                        if len(record.seq)<max(gf[lastex])+3:
                            reachtheend=True
                            break
                        newcodon=str(record.seq[max(gf[lastex]):max(gf[lastex])+3])
                        gf[lastex]={min(gf[lastex]),max(gf[lastex])+3}
                        nlrsq=nlrsq+newcodon
                        prot=prot+str(Seq(newcodon).translate())
                    elif hits[g]["sens"]=="-":
                        if min(gf[lastex])-4<0:
                            reachtheend=True
                            break
                        newcodon=str(record.seq[min(gf[lastex])-4:min(gf[lastex])-1].reverse_complement())
                        gf[lastex]={min(gf[lastex])-3,max(gf[lastex])}
                        nlrsq=nlrsq+newcodon
                        prot=prot+str(Seq(newcodon).translate())
    if prot[0]!="M":
        endloop=False
        while endloop==False:
            for record in SeqIO.parse(args.assembly_file,"fasta"):
                if str(record.id)==hits[g]["scaffold"]:
                    if hits[g]["sens"]=="+":
                        if min(gf[0])-3<1:
                            endloop=True
                            break
                        newcodon=str(record.seq[min(gf[0])-4:min(gf[0])-1])
                        if str(Seq(newcodon).translate())=="*":
                            endloop=True
                            break
                        else:
                            gf[0]={min(gf[0])-3,max(gf[0])}
                            nlrsq=newcodon+nlrsq
                            prot=str(Seq(newcodon).translate())+prot
                            if str(Seq(newcodon).translate())=="M":
                                endloop=True
                                break
                    elif hits[g]["sens"]=="-":
                        if max(gf[0])+3>len(record.seq):
                            endloop=True
                            break
                        newcodon=str(record.seq[max(gf[0]):max(gf[0])+3])
                        if str(Seq(newcodon).translate().reverse_complement())=="*":
                            endloop=True
                            break
                        else:
                            gf[0]={min(gf[0]),max(gf[0])+3}
                            nlrsq=newcodon+nlrsq
                            prot=str(Seq(newcodon).translate().reverse_complement())+prot
                            if str(Seq(newcodon).translate().reverse_complement())=="M":
                                endloop=True
                                break
        if prot[0]!="M":
            if prot.find("M")==-1:
                if g in nlrlist_prolong:
                    nlrlist_prolong.remove(g)
                elif g in sacerlist_kept:
                    sacerlist_kept.remove(g)
                elif g in newgenelist:
                    newgenelist.remove(g)
                continue
            else:
                toremove=prot.find("M")*3
                nlrsq=nlrsq[toremove:]
                prot=prot[int(toremove/3):]
                if hits[g]["sens"]=="+":
                    indices=sorted(list(gf.keys())).copy()
                    for i in indices:
                        if toremove>=(max(gf[i])-min(gf[i])+1):
                            toremove=toremove-(max(gf[i])-min(gf[i])+1)
                            gf.pop(i)
                            if toremove==0:
                                gf=sortgff(gf,hits[g]["sens"])
                                break
                        else:
                            gf[i]={min(gf[i])+toremove,max(gf[i])}
                            break
                elif hits[g]["sens"]=="-":
                    indices=sorted(list(gf.keys())).copy()
                    for i in indices:
                        if toremove>=(max(gf[i])-min(gf[i])+1):
                            toremove=toremove-(max(gf[i])-min(gf[i])+1)
                            gf.pop(i)
                            if toremove==0:
                                gf=sortgff(gf,hits[g]["sens"])
                                break
                        else:
                            gf[i]={min(gf[i]),max(gf[i])-toremove}
                            break
    l=len(nlrsq)
    if l<int(args.minsize_nlr):
        if g in nlrlist_prolong:
            nlrlist_prolong.remove(g)
        elif g in sacerlist_kept:
            sacerlist_kept.remove(g)
        elif g in newgenelist:
            newgenelist.remove(g)
    else:
        hits[g]["gff"]=gf
        fastant_hits[g]=nlrsq
        fastaaa_hits[g]=prot

# once again we do the same, but this time for "unknown" hits
# a difference is that here, we also store the gffs in a new dico, named gff_unk
# because here we name the unknown genes as ind+"_newnlr"+str(k)
fastant_unk={}
fastaaa_unk={}
gff_unk={}
score_unk={}
sens_unk={}
min_unk={}
max_unk={}
k=0
for scaffold in hits["unknown"]:
    keylist=set(hits["unknown"][scaffold].keys()).copy()
    for j in keylist:
        nlrsq=""
        gf=sortgff(hits["unknown"][scaffold][j]["gff"],hits["unknown"][scaffold][j]["sens"])
        for i in range(len(gf)):
            for record in SeqIO.parse(args.assembly_file,"fasta"):
                if str(record.id)==scaffold:
                    if len(gf[i])>1:
                        sq=record.seq[min(gf[i])-1:max(gf[i])]
                    else:
                        sq=Seq(record.seq[min(gf[i])-1])
                    if hits["unknown"][scaffold][j]["sens"]=="-":
                        sq=sq.reverse_complement()
                    nlrsq=nlrsq+str(sq)
        prot=str(Seq(nlrsq).translate())
        if "*" in prot[:-1]:
            print("problem : there is a stop codon in the sequence of unknown hit "+ind+"_newnlr"+str(i))
            print("prot sequence : "+prot)
            print("nt sequence : "+nlrsq)
            print("gff : ")
            print(gf)
            os.kill()
        if len(nlrsq)%3!=0:
            lastex=max(gf.keys())
            missing=3-(len(nlrsq)%3)
            for record in SeqIO.parse(args.assembly_file,"fasta"):
                if str(record.id)==scaffold:
                    if hits["unknown"][scaffold][j]["sens"]=="+":
                        if len(record.seq)<max(gf[lastex])+3:
                            break
                        newcodon=str(record.seq[max(gf[lastex]):max(gf[lastex])+missing])
                        gf[lastex]={min(gf[lastex]),max(gf[lastex])+missing}
                        nlrsq=nlrsq+newcodon
                        prot=prot+str(Seq(newcodon).translate())
                    elif hits["unknown"][scaffold][j]["sens"]=="-":
                        if min(gf[lastex])-4<0:
                            break
                        newcodon=str(record.seq[min(gf[lastex])-1-missing:min(gf[lastex])-1].reverse_complement())
                        gf[lastex]={min(gf[lastex])-missing,max(gf[lastex])}
                        nlrsq=nlrsq+newcodon
                        prot=prot+str(Seq(newcodon).translate())
        if prot[-1]!="*":
            lastex=max(gf.keys())
            newcodon=""
            reachtheend=False
            while str(Seq(newcodon).translate())!="*" and reachtheend==False:
                for record in SeqIO.parse(args.assembly_file,"fasta"):
                    if str(record.id)==scaffold:
                        if hits["unknown"][scaffold][j]["sens"]=="+":
                            if len(record.seq)<max(gf[lastex])+3:
                                reachtheend=True
                                break
                            newcodon=str(record.seq[max(gf[lastex]):max(gf[lastex])+3])
                            gf[lastex]={min(gf[lastex]),max(gf[lastex])+3}
                            nlrsq=nlrsq+newcodon
                            prot=prot+str(Seq(newcodon).translate())
                        elif hits["unknown"][scaffold][j]["sens"]=="-":
                            if min(gf[lastex])-4<0:
                                reachtheend=True
                                break
                            newcodon=str(record.seq[min(gf[lastex])-4:min(gf[lastex])-1].reverse_complement())
                            gf[lastex]={min(gf[lastex])-3,max(gf[lastex])}
                            nlrsq=nlrsq+newcodon
                            prot=prot+str(Seq(newcodon).translate())
        if prot[0]!="M":
            endloop=False
            while endloop==False:
                for record in SeqIO.parse(args.assembly_file,"fasta"):
                    if str(record.id)==scaffold:
                        if hits["unknown"][scaffold][j]["sens"]=="+":
                            if min(gf[0])-3<1:
                                endloop=True
                                break
                            newcodon=str(record.seq[min(gf[0])-4:min(gf[0])-1])
                            if str(Seq(newcodon).translate())=="*":
                                endloop=True
                                break
                            else:
                                gf[0]={min(gf[0])-3,max(gf[0])}
                                nlrsq=newcodon+nlrsq
                                prot=str(Seq(newcodon).translate())+prot
                                if str(Seq(newcodon).translate())=="M":
                                    endloop=True
                                    break
                        elif hits["unknown"][scaffold][j]["sens"]=="-":
                            if max(gf[0])+3>len(record.seq):
                                endloop=True
                                break
                            newcodon=str(record.seq[max(gf[0]):max(gf[0])+3])
                            if str(Seq(newcodon).translate().reverse_complement())=="*":
                                endloop=True
                                break
                            else:
                                gf[0]={min(gf[0]),max(gf[0])+3}
                                nlrsq=newcodon+nlrsq
                                prot=str(Seq(newcodon).translate().reverse_complement())+prot
                                if str(Seq(newcodon).translate().reverse_complement())=="M":
                                    endloop=True
                                    break
            if prot[0]!="M":
                if prot.find("M")==-1:
                    hits["unknown"][scaffold].pop(j)
                    continue
                else:
                    toremove=prot.find("M")*3
                    nlrsq=nlrsq[toremove:]
                    prot=prot[int(toremove/3):]
                    if hits["unknown"][scaffold][j]["sens"]=="+":
                        indices=sorted(list(gf.keys())).copy()
                        for i in indices:
                            if toremove>=(max(gf[i])-min(gf[i])+1):
                                toremove=toremove-(max(gf[i])-min(gf[i])+1)
                                gf.pop(i)
                                if toremove==0:
                                    gf=sortgff(gf,hits["unknown"][scaffold][j]["sens"])
                                    break
                            else:
                                gf[i]={min(gf[i])+toremove,max(gf[i])}
                                break
                    elif hits["unknown"][scaffold][j]["sens"]=="-":
                        indices=sorted(list(gf.keys())).copy()
                        for i in indices:
                            if toremove>=(max(gf[i])-min(gf[i])+1):
                                toremove=toremove-(max(gf[i])-min(gf[i])+1)
                                gf.pop(i)
                                if toremove==0:
                                    gf=sortgff(gf,hits["unknown"][scaffold][j]["sens"])
                                    break
                            else:
                                gf[i]={min(gf[i]),max(gf[i])-toremove}
                                break
        l=len(nlrsq)
        if l<int(args.minsize_nlr):
            hits["unknown"][scaffold].pop(j)
        else:
            hits["unknown"][scaffold][j]["gff"]=gf
            fastant_unk[ind+"_newnlr"+str(k)]=nlrsq
            fastaaa_unk[ind+"_newnlr"+str(k)]=prot
            if scaffold not in gff_unk:
                gff_unk[scaffold]={}
            gff_unk[scaffold][ind+"_newnlr"+str(k)]=gf
            score_unk[ind+"_newnlr"+str(k)]=hits["unknown"][scaffold][j]["bestscore"]
            sens_unk[ind+"_newnlr"+str(k)]=hits["unknown"][scaffold][j]["sens"]
            min_unk[ind+"_newnlr"+str(k)]=hits["unknown"][scaffold][j]["debut"]
            max_unk[ind+"_newnlr"+str(k)]=hits["unknown"][scaffold][j]["fin"]
            k+=1
        
with open("./readable_"+ind+"_correc_exo.txt",'a') as readable: 
    readable.write("#sequences extracted and checked\n")
    readable.write("#time is "+str(timedelta(seconds=(timer()-t0)))+"\n#\n")
    readable.write("#\n#calcul of statistics\n")
  
###############################################################################
######## 14. writing of the outputs
###############################################################################

with open("./readable_"+ind+"_correc_exo.txt",'a') as readable:    
    readable.write("# "+str(len(fastant_nlrnotouch))+" NLRs are not touched (i.e. not prolongated)\n")
    readable.write(str(set(nlrlist_notouch.keys()))+"\n")
    readable.write("#\n")
    readable.write("# we identified "+str(len(nlrlist_prolong)+len(sacerlist_kept)+len(newgenelist))+
                   " previous genes and "+
          str(sum([len(hits["unknown"][x]) for x in hits["unknown"]]))+" new genes as NLRs\n")
    readable.write("# among them are "+str(len(nlrlist_prolong))+" prolongated NLRs, "+
                   str(len(sacerlist_kept))+
          " sacer genes, and "+str(len(newgenelist))+" genesnot previously characterized as NLR\n")
    readable.write("# "+str(len(sacer_notkept))+" sacer genes were not kept, or only as pseudogenes\n")
    readable.write("# "+str(len(nlrnonb))+" NLR genes did not have a NB signal on a total of "+
                   str(nlr_ind)+" NLR genes, so they we removed\n")
    readable.write("# "+str(len(nlrlist_pseudo))+" other NLR genes did have a bad NB signal or"+
                   " an incorrect structure, so they are considered as pseudogenes\n")
    readable.write("#\n")
    readable.write("# all the former NLRs of "+ind+"\n")
    readable.write(str(nlrlist_ind)+"\n")
    readable.write("#\n")
    readable.write("# NLRs that are prolongated\n")
    readable.write(str(nlrlist_prolong)+"\n")
    readable.write("#\n")
    readable.write("# NLRs that are pseudogenes\n")
    readable.write(str(nlrlist_pseudo)+"\n")
    readable.write("#\n")
    readable.write("# NLRs that are removed\n")
    readable.write(str(nlrnonb)+"\n")
    readable.write("#\n")
    readable.write("# sacers identified as NLRs\n")
    readable.write(str(sacerlist_kept)+"\n")
    readable.write("#\n")
    readable.write("# new genes identified as NLRs\n")
    readable.write(str(newgenelist)+"\n")
    readable.write("#\n")
    readable.write("#\n")

with open("./readable_"+ind+"_correc_exo.txt",'a') as readable:
    
    if args.output is not None:
        with open(args.output,"w") as OUT:
            OUT.write("#analyse exonerate of "+ind+" based on a set of "+str(len(nlrs))+" NLR genes\n")
            OUT.write("# NLRs not prolongated : "+" ".join(list(nlrlist_notouch.keys()))+"\n#\n")
            OUT.write("#gene\tscore_exonerate\tscaffold\tstrand\tstart\tend\tgff\n")
            OUT.write("# prolongated NLRs\n")
            for g in nlrlist_prolong:
                OUT.write(g+"\t"+str(hits[g]["bestscore"])+"\t"+str(hits[g]["scaffold"])+"\t"+
                          str(hits[g]["sens"])+"\t"+
                          str(hits[g]["debut"])+"\t"+str(hits[g]["fin"])+"\t"+str(hits[g]["gff"])+"\n")
            OUT.write("# sacers kept as NLR\n")
            for g in sacerlist_kept:
                OUT.write(g+"\t"+str(hits[g]["bestscore"])+"\t"+str(hits[g]["scaffold"])+"\t"+
                          str(hits[g]["sens"])+"\t"+
                          str(hits[g]["debut"])+"\t"+str(hits[g]["fin"])+"\t"+str(hits[g]["gff"])+"\n")
            OUT.write("# new genes identified as NLRs\n")
            for g in newgenelist:
                OUT.write(g+"\t"+str(hits[g]["bestscore"])+"\t"+str(hits[g]["scaffold"])+"\t"+
                          str(hits[g]["sens"])+"\t"+
                          str(hits[g]["debut"])+"\t"+str(hits[g]["fin"])+"\t"+str(hits[g]["gff"])+"\n")
            OUT.write("# unknown genes identified as NLRs\n")
            for scaffold in gff_unk:
                for g in gff_unk[scaffold].keys():
                    OUT.write(g+"\t"+str(score_unk[g])+
                              "\t"+scaffold+"\t"+str(sens_unk[g])+"\t"+
                              str(min_unk[g])+"\t"+
                              str(max_unk[g])+"\t"+
                              str(gff_unk[scaffold][g])+"\n")
            OUT.write("# there was also "+str(nb_pseudo)+" pseudo-NLRs identified\n")
            OUT.write("# "+str(len(nlrlist_pseudo))+" NLRs are now pseudo-NLRs\n")
            OUT.write("# "+str(nlrlist_pseudo)+"\n")
            OUT.write("# all the pseudo-NLRs identified\n")
            j=0
            for scaffold in sorted(list(pseudonlr.keys())):
                for i in sorted(list(pseudonlr[scaffold].keys())):
                    OUT.write(ind+"_"+"pseudo"+str(j)+"\t"+str(pseudonlr[scaffold][i]["score"])+
                              "\t"+scaffold+"\t"+
                              str(pseudonlr[scaffold][i]["debut"])+"\t"+
                              str(pseudonlr[scaffold][i]["fin"])+"\t"+
                              str(pseudonlr[scaffold][i]["genes"])+"\n")
    
    if args.fasta_nlrs_nt_out is not None:
        with open(args.fasta_nlrs_nt_out,"w") as OUT:
            for g in fastant_nlrnotouch:
                OUT.write(">"+g+"\n"+fastant_nlrnotouch[g]+"\n")
            for g in nlrlist_prolong:
                OUT.write(">"+g+"\n"+fastant_hits[g]+"\n")
            for g in sacerlist_kept:
                OUT.write(">"+g+"\n"+fastant_hits[g]+"\n")
            for g in newgenelist:
                OUT.write(">"+g+"\n"+fastant_hits[g]+"\n")
            for g in fastant_unk:
                OUT.write(">"+g+"\n"+fastant_unk[g]+"\n")
    if args.fasta_nlrs_aa_out is not None:
        with open(args.fasta_nlrs_aa_out,"w") as OUT:
            for g in fastaaa_nlrnotouch:
                OUT.write(">"+g+"\n"+fastaaa_nlrnotouch[g]+"\n")
            for g in nlrlist_prolong:
                OUT.write(">"+g+"\n"+fastaaa_hits[g]+"\n")
            for g in sacerlist_kept:
                OUT.write(">"+g+"\n"+fastaaa_hits[g]+"\n")
            for g in newgenelist:
                OUT.write(">"+g+"\n"+fastaaa_hits[g]+"\n")
            for g in fastant_unk:
                OUT.write(">"+g+"\n"+fastaaa_unk[g]+"\n")
    if args.pseudo_out is not None:
        with open(args.pseudo_out,"w") as OUT:
            j=0
            for scaffold in sorted(list(pseudonlr.keys())):
                for i in sorted(list(pseudonlr[scaffold].keys())):
                    OUT.write(ind+"_"+"pseudo"+str(j)+"\t"+str(pseudonlr[scaffold][i]["score"])+
                              "\t"+scaffold+"\t"+
                              str(pseudonlr[scaffold][i]["debut"])+"\t"+
                              str(pseudonlr[scaffold][i]["fin"])+"\t"+
                              str(pseudonlr[scaffold][i]["genes"])+"\n")
    
    readable.write("#complete running time is "+str(timedelta(seconds=(timer()-t0)))+"\n#\n")
            
















