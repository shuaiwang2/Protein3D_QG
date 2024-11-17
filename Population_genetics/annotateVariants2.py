#!python

import NucleotideCodeSubstitution
import FastaFile
import GffFile
import MyUtil
import math
import copy
import sys
# bs674@cornell.edu

_buckets = []

def read_data(gffFile, fastaFile, vcfFilePath):
    all_canonicalTranscripts = set()
    with open("/media/shuaiwang/3d/NAMpopulation_alphafold/allNAMs/canonicalPorteins/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.canonical_transcripts") as f:
        for line in f:
            element = line.split()
            all_canonicalTranscripts.add(element[0])

    chromosome_transcript_dict, chromosome_transcript_list = GffFile.readGff( gffFile, all_canonicalTranscripts )
    chromosome_names, fastas = FastaFile.readFastaFile( fastaFile )
    GffFile.update_sequence_information(fastas, chromosome_transcript_dict)
    with open(vcfFilePath) as f:
        next(f)
        for line in f:
            element = line.split()
            chr = element[0]
            position = int(element[1])
            ref = element[4].split(':')
            alt = element[5].split(':')
            ref_allele = ref[0]
            ref_frq = float(ref[1])
            alt_allele = alt[0]
            alt_frq = float(alt[1])
#            print(ref_allele + "\t" + fastas[chr][position-1] + "\t" + alt_allele)
            if alt_frq >0 and ref_frq >0 and (chr in {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"}):
                assert ref_allele == fastas[chr][position-1], "not working"
                transcript_name = MyUtil.overlap_with_certain_transcript(position, chr, chromosome_transcript_dict)
                if None != transcript_name:
                    if "+" == chromosome_transcript_dict[chr][transcript_name].strand:
                        thisCdsPosition = 1
                        for cds in chromosome_transcript_dict[chr][transcript_name].Cds:
                            if cds[0]<=position and position <= cds[1]:
                                thisCdsPosition += (position-cds[0])
                                amino_acid_index = math.ceil((thisCdsPosition)/3.0)
                                alterCdsSequence = chromosome_transcript_dict[chr][transcript_name].cds_sequence[0:thisCdsPosition-1] + alt_allele + chromosome_transcript_dict[chr][transcript_name].cds_sequence[thisCdsPosition:]
                                codon = chromosome_transcript_dict[chr][transcript_name].cds_sequence[(amino_acid_index-1)*3:((amino_acid_index-1)*3+3)]
                                altCodon = alterCdsSequence[(amino_acid_index-1)*3:((amino_acid_index-1)*3+3)]
                                if NucleotideCodeSubstitution.middleStandardGeneticCode[codon] != NucleotideCodeSubstitution.middleStandardGeneticCode[altCodon]:
                                    print("nonsynonymous" + "\t" + chr + "\t" + str(position) + "\t" + element[3] + "\t" + transcript_name + "\t" + ref_allele + "\t" + alt_allele + "\t" + str(amino_acid_index) + "\t" + str(ref_frq) + "\t" + str(alt_frq) + "\t" + codon + "\t" + altCodon + "\t" + str(thisCdsPosition))
                                else:
                                    print("synonymous" + "\t" + chr + "\t" + str(position) + "\t" + element[3] + "\t" + transcript_name + "\t" + ref_allele + "\t" + alt_allele + "\t" + str(amino_acid_index) + "\t" + str(ref_frq) + "\t" + str(alt_frq) + "\t" + codon + "\t" + altCodon + "\t" + str(thisCdsPosition))
                            else:
                                thisCdsPosition = thisCdsPosition + cds[1] - cds[0] + 1
                    else:
                        thisCdsPosition = 1
                        for cds in chromosome_transcript_dict[chr][transcript_name].Cds:
                            if cds[0]<=position and position <= cds[1]:
                                thisCdsPosition += (cds[1]-position)
                                amino_acid_index = math.ceil((thisCdsPosition)/3.0)
                                alterCdsSequence = chromosome_transcript_dict[chr][transcript_name].cds_sequence[0:thisCdsPosition-1] + NucleotideCodeSubstitution.getComplementary(alt_allele) + chromosome_transcript_dict[chr][transcript_name].cds_sequence[thisCdsPosition:]
                                codon = chromosome_transcript_dict[chr][transcript_name].cds_sequence[(amino_acid_index-1)*3:((amino_acid_index-1)*3+3)]
                                altCodon = alterCdsSequence[(amino_acid_index-1)*3:((amino_acid_index-1)*3+3)]
                                if NucleotideCodeSubstitution.middleStandardGeneticCode[codon] != NucleotideCodeSubstitution.middleStandardGeneticCode[altCodon]:
                                    print("nonsynonymous" + "\t"+ chr + "\t" + str(position) + "\t" + element[3] + "\t" + transcript_name + "\t" + ref_allele + "\t" + alt_allele + "\t" + str(amino_acid_index) + "\t" + str(ref_frq) + "\t" + str(alt_frq) + "\t" + codon + "\t" + altCodon + "\t" + str(thisCdsPosition))
                                else:
                                    print("synonymous" + "\t" + chr + "\t" + str(position) + "\t" + element[3] + "\t" + transcript_name + "\t" + ref_allele + "\t" + alt_allele + "\t" + str(amino_acid_index) + "\t" + str(ref_frq) + "\t" + str(alt_frq) + "\t" + codon + "\t" + altCodon + "\t" + str(thisCdsPosition))

                            else:
                                thisCdsPosition = thisCdsPosition + cds[1] - cds[0] + 1

    return

read_data("/media/shuaiwang/3d/282/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3", "/media/shuaiwang/3d/282vcf_phen_song/282_libs_2015/uplifted_APGv4/Zm-B73-REFERENCE-NAM-5.0.fa", sys.argv[1])
