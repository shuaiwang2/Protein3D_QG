import numpy as np
import re
import sys
import pickle
import glob
import csv
import os
# bs674@cornell.edu

# this function was copied from the alphafold source code

def readMsa(fastaFile, uniqueSeq = True):
    if uniqueSeq:
        seqToName = {}
        name = ""
        seq = []
        with open(fastaFile) as f:
            for line in f:
                m = re.search('^>(\S+)', line)
                if (m != None):
                    if (len(name) > 0) & (len(seq) > 0):
                        s = ''.join(seq)
                        s = re.sub("\\s", "", s)
                        s = s.upper()
                        if s in seqToName:
                            seqToName[s] = seqToName[s] + "," + name
                        else:
                            seqToName[s] = name
                    name = m.group(1)
                    seq = []
                else:
                    seq.append(line)
            if (len(name) > 0) & (len(seq) > 0):
                s = ''.join(seq)
                s = re.sub("\\s", "", s)
                s = s.upper()
                if s in seqToName:
                    seqToName[s] = seqToName[s] + "," + name
                else:
                    seqToName[s] = name
        seq_names = []
        fastas = {}
        for seq in sorted(seqToName):
            seq_names.append(seqToName[seq])
            fastas[seqToName[seq]]=seq
        return seq_names, fastas
    else:
        fastas = {}
        seq_names = []
        name = ""
        seq = []
        with open(fastaFile) as f:
            for line in f:
                m = re.search('^>(\S+)', line)
                if (m != None):
                    if (len(name) > 0) & (len(seq) > 0):
                        s = ''.join(seq)
                        s = re.sub("\\s", "", s)
                        s = s.upper()
                        fastas[name] = s
                        seq_names.append(name)
                    name = m.group(1)
                    seq = []
                else:
                    seq.append(line)
            if (len(name) > 0) & (len(seq) > 0):
                s = ''.join(seq)
                s = re.sub("\\s", "", s)
                s = s.upper()
                fastas[name] = s
                seq_names.append(name)
        return seq_names, fastas

class Coordinate:
    x = 0
    y = 0
    z = 0
    plddtScore = 0
    def __init__(self, x, y, z, plddtScore):
        self.x = x
        self.y = y
        self.z = z
        self.plddtScore = plddtScore

def readPdb(pdbFile):
    coordinates = {}
    with open(pdbFile) as f:
        for line in f:
            elements = line.split()
            if elements[0] == "ATOM":
                atomId = elements[2]
                aaid = elements[5]
                x = elements[6]
                y = elements[7]
                z = elements[8]
##                print(line)
                plddtScore = float(elements[10])
                coordinate = Coordinate(x, y, z, plddtScore)
                if aaid not in coordinates:
                    coordinates[aaid] = {}
                coordinates[aaid][atomId] = coordinate
    return coordinates

def calculateDistance(pdbFileRef, pdbFileQuery, positions):
    coordinatesRef = readPdb(pdbFileRef)
    coordinatesQuery = readPdb(pdbFileQuery)
    r = []
    q = []
    for refPosition in sorted(positions):
        if coordinatesRef[str(refPosition)]["CA"].plddtScore >= 70.0 and coordinatesQuery[str(positions[refPosition])]["CA"].plddtScore >= 70.0:
            r.append([coordinatesRef[str(refPosition)]["CA"].x, coordinatesRef[str(refPosition)]["CA"].y, coordinatesRef[str(refPosition)]["CA"].z])
            q.append([coordinatesQuery[str(positions[refPosition])]["CA"].x, coordinatesQuery[str(positions[refPosition])]["CA"].y, coordinatesQuery[str(positions[refPosition])]["CA"].z])

            # r.append([coordinatesRef[str(refPosition)]["N"].x, coordinatesRef[str(refPosition)]["N"].y, coordinatesRef[str(refPosition)]["N"].z])
            # q.append([coordinatesQuery[str(positions[refPosition])]["N"].x, coordinatesQuery[str(positions[refPosition])]["N"].y, coordinatesQuery[str(positions[refPosition])]["N"].z])
            #
            # r.append([coordinatesRef[str(refPosition)]["C"].x, coordinatesRef[str(refPosition)]["C"].y, coordinatesRef[str(refPosition)]["C"].z])
            # q.append([coordinatesQuery[str(positions[refPosition])]["C"].x, coordinatesQuery[str(positions[refPosition])]["C"].y, coordinatesQuery[str(positions[refPosition])]["C"].z])
            #
            # r.append([coordinatesRef[str(refPosition)]["O"].x, coordinatesRef[str(refPosition)]["O"].y, coordinatesRef[str(refPosition)]["O"].z])
            # q.append([coordinatesQuery[str(positions[refPosition])]["O"].x, coordinatesQuery[str(positions[refPosition])]["O"].y, coordinatesQuery[str(positions[refPosition])]["O"].z])
    return lddt(np.array([q], dtype=np.float32), np.array([r], dtype=np.float32), np.array([[[1]] * len(r)], dtype=np.float32), cutoff=15., per_residue=False)

maxInt = sys.maxsize
csv.field_size_limit(maxInt)
def readRilPanGene(rilPanGeneFile, pan_gene_id):
    ril_panGenesId_dict = dict()
    ril_accession_dict = dict()
    f1 = open(rilPanGeneFile, 'r')
    i = 0
    for line in f1:
        line = line.strip()
        l = line.split()
        if i == 0:
            ii = -1
            for l1 in l:
                if ii >= 0:
                    ril_accession_dict[l1] = ii
                ii = ii + 1
        elif i > 0:
            pan_gene = l[0]
            ril_panGenesId_dict[pan_gene] = i - 1
        i = i+1
    f1.close()
    row_index = ril_panGenesId_dict[pan_gene_id] + 1
    with open(rilPanGeneFile, 'r') as fin:
        reader=csv.reader(fin)
        result=[[s for s in row] for i,row in enumerate(reader) if i == row_index]
    return ril_accession_dict, result[0][0].split()
def blosum62Score():
    aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']
    scores = {}
    blosumMatrix = [[4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4],
                    [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4],
                    [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4],
                    [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4],
                    [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4],
                    [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4],
                    [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4],
                    [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4],
                    [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4],
                    [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4],
                    [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4],
                    [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4],
                    [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4],
                    [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4],
                    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4],
                    [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4],
                    [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4],
                    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4],
                    [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4],
                    [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4],
                    [-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4],
                    [-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4],
                    [0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4],
                    [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1]]

    for i in range(len(aa)):
        scores[aa[i]] = {}
        for j in range(len(aa)):
            scores[aa[i]][aa[j]] = blosumMatrix[i][j]
    return scores

def getDistance(a, b, scores):
    a = re.sub("\\s", "", a)
    b = re.sub("\\s", "", b)
    a = a.upper()
    b = b.upper()
    positives = 0
    for i in range(len(a)):
        if (not (a[i]=='-' or b[i]=='-')) and scores[a[i]][b[i]] > 0:
            positives = positives + 1
    a = re.sub("-", "", a)
    b = re.sub("-", "", b)
    return (positives*2)/(len(a) + len(b))

if __name__ == '__main__':
    fastas = pickle.load(open("/media/bs674/ppi8t/NAMpopulation_alphafold/allNAMs/canonicalPorteins/shapmer/fastas", "rb"))
    fastasRev = pickle.load(open("/media/bs674/ppi8t/NAMpopulation_alphafold/allNAMs/canonicalPorteins/shapmer/fastasRev", "rb"))
    fastasRev2 = pickle.load(open("/media/bs674/ppi8t/NAMpopulation_alphafold/allNAMs/canonicalPorteins/shapmer/fastasRev2", "rb"))

    proteinId2Accession = dict()
    f1 = open("/media/bs674/ppi8t/NAMpopulation_alphafold/allNAMs/ProteinId_to_accession", 'r')
    for line in f1:
        line = line.strip()
        l = line.split()
        proteinId2Accession[l[0]] = l[1]

    pan_group_id = sys.argv[1]

    seq_names, msaFastas = readMsa("/media/bs674/ppi8t/NAMpopulation_alphafold/allNAMs/canonicalPorteins/deleteriousMutations/" + pan_group_id + "/sequences.mafft", uniqueSeq = False)
    similarity_matrix = np.full(shape=(len(seq_names), len(seq_names)), fill_value=1.0, dtype=float)
    scores = blosum62Score()
    for i, seq_name_i in enumerate(seq_names):
        for j, seq_name_j in enumerate(seq_names):
            if i < j:
                assert len(msaFastas[seq_name_i]) == len(msaFastas[seq_name_j]), "sequence length error"
                similarity_matrix[i][j] = getDistance(msaFastas[seq_name_i], msaFastas[seq_name_j], scores)
                similarity_matrix[j][i] = similarity_matrix[i][j]

    output = open("/media/bs674/ppi8t/NAMpopulation_alphafold/allNAMs/canonicalPorteins/deleteriousMutations/" + pan_group_id + "/mafftSequenceSimilarityMatrix", 'w')

    for row_label, row in zip(seq_names, similarity_matrix):
        output.write('%s %s' % (row_label, ' '.join('%03s' % i for i in row)) + "\n")
        #output.write ( '%s %s' % (row_label, ' '.join('%03s' % i for i in row)) + "\n" )
    output.close()
