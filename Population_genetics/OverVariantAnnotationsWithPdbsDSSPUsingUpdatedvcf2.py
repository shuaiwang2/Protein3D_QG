#!python

import NucleotideCodeSubstitution
import FastaFile
import re
import GffFile
import MyUtil
import math
import copy
import sys
import glob
#from Bio.SubsMat import MatrixInfo
from Bio.Align import substitution_matrices



# bs674@cornell.edu

def score_match(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]
_buckets = []

blosum = substitution_matrices.load('BLOSUM62') 
#blosum = MatrixInfo.blosum62

def read_data():
    protein_names, fastas = FastaFile.readFastaFile( "/media/shuaiwang/3d/NAMpopulation_alphafold/allNAMs/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa" )
    seq_to_id={}
    name = ""
    seq = []
    with open("/media/shuaiwang/3d/onedrive/OneDrive-2024-03-13/B73V5/unique.seqs.fa") as f:
        for line in f:
            m = re.search('^>(\S+)', line)
            if (m != None):
                if (len(name) > 0) & (len(seq) > 0):
                    s = ''.join(seq)
                    s = re.sub("\\s", "", s)
                    s = s.upper()
                    seq_to_id[s]=name
                name = m.group(1)
                seq = []
            else:
                seq.append(line)
        if (len(name) > 0) & (len(seq) > 0):
            s = ''.join(seq)
            s = re.sub("\\s", "", s)
            s = s.upper()
            seq_to_id[s]=name

    wanted_dict={}
    other_zea_variants = {}
    with open("/media/shuaiwang/3d/282/all.freq") as f:
        for line in f:
            e = line.split()
            if e[1] not in wanted_dict:
                wanted_dict[e[1]] = set()
                other_zea_variants[e[1]] = {}
            wanted_dict[e[1]].add(e[2])
            other_zea_variants[e[1]][int(e[2])] = ["NA", "NA","NA", "NA","NA"]

    with open("/media/shuaiwang/3d/outgroup/outgroupAlignment/outgroup/TIL11.gvcf") as f:
        for line in f:
            if line[0] != '#':
                element = line.split()
                start = int(element[1])
                if element[4][0] == '<':
                    end = int(element[7].split(";END=")[1])
                    for i in range(start, end+1):
                        if (element[0] in other_zea_variants) and (i in other_zea_variants[element[0]]):
                            other_zea_variants[element[0]][i][0] = "R"
                            other_zea_variants[element[0]][i][1] = "R"
#                            print(element[0] + "\t" + str(i) + "\t" + other_zea_variants[element[0]][i])
                elif (element[0] in other_zea_variants) and (start in other_zea_variants[element[0]]):
                    other_zea_variants[element[0]][start][0] = element[3]
                    other_zea_variants[element[0]][start][1] = element[4].split(",")[0]

    with open("/media/shuaiwang/3d/outgroup/outgroupAlignment/outgroup/TIL18.gvcf") as f:
        for line in f:
            if line[0] != '#':
                element = line.split()
                start = int(element[1])
                if element[4][0] == '<':
                    end = int(element[7].split(";END=")[1])
                    for i in range(start, end+1):
                        if (element[0] in other_zea_variants) and (i in other_zea_variants[element[0]]):
                            other_zea_variants[element[0]][i][2] = "R"
                            if other_zea_variants[element[0]][i][0] == "NA":
                                other_zea_variants[element[0]][i][0] = "R"
                elif (element[0] in other_zea_variants) and (start in other_zea_variants[element[0]]):
                    other_zea_variants[element[0]][start][2] = element[4].split(",")[0]
                    if other_zea_variants[element[0]][start][0] == "NA":
                        other_zea_variants[element[0]][start][0] = element[3]

    with open("/media/shuaiwang/3d/outgroup/outgroupAlignment/outgroup/luxurians.gvcf") as f:
        for line in f:
            if line[0] != '#':
                element = line.split()
                start = int(element[1])
                if element[4][0] == '<':
                    end = int(element[7].split(";END=")[1])
                    for i in range(start, end+1):
                        if (element[0] in other_zea_variants) and (i in other_zea_variants[element[0]]):
                            other_zea_variants[element[0]][i][3] = "R"
                            if other_zea_variants[element[0]][i][0] == "NA":
                                other_zea_variants[element[0]][i][0] = "R"
                elif (element[0] in other_zea_variants) and (start in other_zea_variants[element[0]]):
                    other_zea_variants[element[0]][start][3] = element[4].split(",")[0]
                    if other_zea_variants[element[0]][start][0] == "NA":
                        other_zea_variants[element[0]][start][0] = element[3]

    with open("/media/shuaiwang/3d/outgroup/outgroupAlignment/outgroup/diploperennis.gvcf") as f:
        for line in f:
            if line[0] != '#':
                element = line.split()
                start = int(element[1])
                if element[4][0] == '<':
                    end = int(element[7].split(";END=")[1])
                    for i in range(start, end+1):
                        if (element[0] in other_zea_variants) and (i in other_zea_variants[element[0]]):
                            other_zea_variants[element[0]][i][4] = "R"
                            if other_zea_variants[element[0]][i][0] == "NA":
                                other_zea_variants[element[0]][i][0] = "R"
                elif (element[0] in other_zea_variants) and (start in other_zea_variants[element[0]]):
                    other_zea_variants[element[0]][start][4] = element[4].split(",")[0]
                    if other_zea_variants[element[0]][start][0] == "NA":
                        other_zea_variants[element[0]][start][0] = element[3]

    gerpdict = {}
    with open("/media/shuaiwang/3d/outgroup/doi_10_5061_dryad_70t85k2__v20181212/GERP/GERP/out_bed.bed") as f:
        for line in f:
            e = line.split()
            if (e[0] in wanted_dict) and (e[1] in wanted_dict[e[0]]):
            #if (e[0] in wanted_dict) and (e[2] in wanted_dict[e[0]]):
                if e[0] not in gerpdict:
                    gerpdict[e[0]] = {}
                gerpdict[e[0]][e[1]] = e[3]

    with open("/media/shuaiwang/3d/282/all.freq") as f:
        for line in f:
   #         print(line)
            [type, chr, position, allele_count, transcript_name, ref_allele, alt_allele, amino_acid_index, ref_frq, alt_frq, codon, altCodon, thisCdsPosition] = line.split()
            protein_name = transcript_name.replace("_T", "_P")
            if (protein_name in protein_names) and (fastas[protein_name] in seq_to_id):
                pdbFile = glob.glob("/media/shuaiwang/3d/NAMpopulation_alphafold/predictedStructures/B73_result/" + seq_to_id[fastas[protein_name]] + "/*mkdssp")
  #              print(pdbFile)
                if len(pdbFile) >0 :
                    pdbFile = pdbFile[0]
                    total_asa = 0.0
                    with open(pdbFile) as pdbf:
                        for pdbline in pdbf:
                            e = pdbline.split()
                            if(len(e) > 10 and e[0] == e[1]) and e[2]=="A" and int(e[0])==int(amino_acid_index):
#                                print(e)
                                total_asa = total_asa + float(pdbline[35:38])
 #                               print(pdbline[13:14].strip() + "\t" + NucleotideCodeSubstitution.middleStandardGeneticCode[codon])
                                assert  pdbline[13:14].strip() == NucleotideCodeSubstitution.middleStandardGeneticCode[codon], "not working"

                    plddt = 0.0
                    pdbFile2 = glob.glob("/media/shuaiwang/3d/NAMpopulation_alphafold/predictedStructures/B73_result/" + seq_to_id[fastas[protein_name]] + "/*maximumplddts.pdb")[0]
                    with open(pdbFile2) as pdbf:
                        for pdbline in pdbf:
                            if (pdbline.startswith('ATOM')) and (int(pdbline[22:26].strip()) == int(amino_acid_index)):
                                plddt = float(pdbline[60:66])
                    #
                    # outAllele = ref_allele
                    # if (chr in other_zea_variants) and (int(position) in other_zea_variants[chr]):
                    #     if other_zea_variants[chr][int(position)] != "R:R":
                    #         outAllele = other_zea_variants[chr][int(position)].split(":")[1]
                    # else:
                    #     outAllele = "NA"

                    if NucleotideCodeSubstitution.middleStandardGeneticCode[codon] != '*':
                        maxAsa = MyUtil.residue_max_acc[MyUtil.singleLetterAminoAcidToTripleLetterAminoAcid[NucleotideCodeSubstitution.middleStandardGeneticCode[codon]]]
                        score = -20
                        if NucleotideCodeSubstitution.middleStandardGeneticCode[altCodon] != '*':
                            pair = (NucleotideCodeSubstitution.middleStandardGeneticCode[codon], NucleotideCodeSubstitution.middleStandardGeneticCode[altCodon])
                            score = score_match(pair, blosum)
                        gerp_score = "NA"
                        if (chr in gerpdict) and (position in gerpdict[chr]):
                            gerp_score = gerpdict[chr][position]
                        print(type + "\t" + chr + "\t" + position + "\t" + allele_count + "\t" + transcript_name + "\t" + seq_to_id[fastas[protein_name]].split("_")[0] + "\t" + seq_to_id[fastas[protein_name]].split("_")[1] + "\t" + ref_allele + "\t" + alt_allele  + "\t" +  other_zea_variants[chr][int(position)][0] + "\t" + other_zea_variants[chr][int(position)][1] + "\t" + other_zea_variants[chr][int(position)][2] + "\t" + other_zea_variants[chr][int(position)][3] + "\t" + other_zea_variants[chr][int(position)][4] + "\t" + amino_acid_index + "\t" + ref_frq + "\t" + alt_frq + "\t" + codon + "\t" +  NucleotideCodeSubstitution.middleStandardGeneticCode[codon] + "\t" + altCodon + "\t" +  NucleotideCodeSubstitution.middleStandardGeneticCode[altCodon] + "\t" + str(score) + "\t" + thisCdsPosition + "\t" + gerp_score + "\t" + str(plddt) + "\t" + str(total_asa) + "\t" + str(maxAsa))
    #                            print(protein_name + "\t" + seq_to_id[fastas[protein_name]] + "\t" + pdbFile )

read_data()
