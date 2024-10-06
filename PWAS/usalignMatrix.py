# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""lDDT protein distance score."""
import numpy as np
import re
import sys
import pickle
import glob
import csv
import os
import re
# bs674@cornell.edu

# this function was copied from the alphafold source code

def lddt(predicted_points,
         true_points,
         true_points_mask,
         cutoff=15.,
         per_residue=False):
    """Measure (approximate) lDDT for a batch of coordinates.

    lDDT reference:
    Mariani, V., Biasini, M., Barbato, A. & Schwede, T. lDDT: A local
    superposition-free score for comparing protein structures and models using
    distance difference tests. Bioinformatics 29, 2722–2728 (2013).

    lDDT is a measure of the difference between the true distance matrix and the
    distance matrix of the predicted points.  The difference is computed only on
    points closer than cutoff *in the true structure*.

    This function does not compute the exact lDDT value that the original paper
    describes because it does not include terms for physical feasibility
    (e.g. bond length violations). Therefore this is only an approximate
    lDDT score.

    Args:
      predicted_points: (batch, length, 3) array of predicted 3D points
      true_points: (batch, length, 3) array of true 3D points
      true_points_mask: (batch, length, 1) binary-valued float array.  This mask
        should be 1 for points that exist in the true points.
      cutoff: Maximum distance for a pair of points to be included
      per_residue: If true, return score for each residue.  Note that the overall
        lDDT is not exactly the mean of the per_residue lDDT's because some
        residues have more contacts than others.

    Returns:
      An (approximate, see above) lDDT score in the range 0-1.
    """

    assert len(predicted_points.shape) == 3
    assert predicted_points.shape[-1] == 3
    assert true_points_mask.shape[-1] == 1
    assert len(true_points_mask.shape) == 3

    # Compute true and predicted distance matrices.
    dmat_true = np.sqrt(1e-10 + np.sum(
        (true_points[:, :, None] - true_points[:, None, :])**2, axis=-1))

    dmat_predicted = np.sqrt(1e-10 + np.sum(
        (predicted_points[:, :, None] -
         predicted_points[:, None, :])**2, axis=-1))

    dists_to_score = (
            (dmat_true < cutoff).astype(np.float32) * true_points_mask *
            np.transpose(true_points_mask, [0, 2, 1]) *
            (1. - np.eye(dmat_true.shape[1]))  # Exclude self-interaction.
    )

    # Shift unscored distances to be far away.
    dist_l1 = np.abs(dmat_true - dmat_predicted)

    # True lDDT uses a number of fixed bins.
    # We ignore the physical plausibility correction to lDDT, though.
    score = 0.25 * ((dist_l1 < 0.5).astype(np.float32) +
                    (dist_l1 < 1.0).astype(np.float32) +
                    (dist_l1 < 2.0).astype(np.float32) +
                    (dist_l1 < 4.0).astype(np.float32))

    # Normalize over the appropriate axes.
    reduce_axes = (-1,) if per_residue else (-2, -1)
    norm = 1. / (1e-10 + np.sum(dists_to_score, axis=reduce_axes))
    score = norm * (1e-10 + np.sum(dists_to_score * score, axis=reduce_axes))
    return score



def readMsa(fastaFile, uniqueSeq = True):
    if uniqueSeq:
        seqToName = {}
        name = ""
        seq = []
        with open(fastaFile) as f:
            for line in f:
                m = re.search('^>(\S+)', line)
                if (m != None) and line[:1] != "#":
                    if (len(name) > 0) & (len(seq) > 0):
                        s = ''.join(seq)
                        s = re.sub("\\s", "", s)
                        s = s.upper()
                        if s in seqToName:
                            seqToName[s] = seqToName[s] + "," + name
                        else:
                            seqToName[s] = name
                    name = m.group(1)
                    name = re.sub(":A", "", name)
                    seq = []
                elif line[:1] != "#":
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
                if (m != None) and line[:1] != "#":
                    if (len(name) > 0) and (len(seq) > 0):
                        s = ''.join(seq)
                        s = re.sub("\\s", "", s)
                        s = s.upper()
                        fastas[name] = s
                        seq_names.append(name)
                    name = m.group(1)
                    name = re.sub(":A", "", name)
                    seq = []
                elif line[:1] != "#":
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
                x = line[30:38]
                y = line[38:46]
                z = line[46:54]
                # print(x)
                # print(y)
                # print(z)
                # print(line)
                plddtScore = float(line[60:66])
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

if __name__ == '__main__':
    pan_group_id = sys.argv[1]
    # seq_names, msaFastas = readMsa("/media/bs674/ppi8t/NAMpopulation_alphafold/allNAMs/canonicalPorteins/deleteriousMutations/" + pan_group_id + "/alignment", uniqueSeq = False)
    seq_names, msaFastas = readMsa("/home/songbs/deleteriousMutations/" + pan_group_id + "/alignment", uniqueSeq = False)
    distance_matrix = np.full(shape=(len(seq_names), len(seq_names)), fill_value=1.0, dtype=float)
    distances = dict()

    distance_matrix_RMSD = np.full(shape=(len(seq_names), len(seq_names)), fill_value=0, dtype=float)
    distances_RMSD = dict()
    for i, seq_name_i in enumerate(seq_names):
        for j, seq_name_j in enumerate(seq_names):
            assert len(msaFastas[seq_name_i]) == len(msaFastas[seq_name_j]), "sequence length error"
            pdbFileRef = "/home/songbs/deleteriousMutations/" + pan_group_id + "/" + seq_name_i
            #pdbFileRef = "/media/bs674/ppi8t/NAMpopulation_alphafold/allNAMs/canonicalPorteins/deleteriousMutations/" + pan_group_id + "/" + seq_name_i
            #pdbFileQuery = "/media/bs674/ppi8t/NAMpopulation_alphafold/allNAMs/canonicalPorteins/deleteriousMutations/" + pan_group_id + "/" + seq_name_j
            pdbFileQuery = "/home/songbs/deleteriousMutations/" + pan_group_id + "/" + seq_name_j
#            print(pdbFileRef)
#            print(pdbFileQuery)
            if pdbFileRef not in distances:
                distances[pdbFileRef] = dict()
                distances_RMSD = dict()

            if pdbFileQuery not in distances:
                distances[pdbFileQuery] = dict()
                distances_RMSD[pdbFileQuery] = dict()

            if pdbFileRef != pdbFileQuery:
                output = os.popen('USalign -mol prot -mm 0 '+ pdbFileRef  + ' ' + pdbFileQuery)
                output = output.read()
                m = re.search('TM-score=\s{0,}(\S+)\s{0,}\(no', output)
                score = m.group(1)
                distance_matrix[i][j] = score
                distances[pdbFileQuery][pdbFileRef] = score
                m1 = re.search('RMSD=\s{0,}(\S+),', output)
                score = m1.group(1)
                distance_matrix_RMSD[i][j] = score
#                distances_RMSD[pdbFileQuery][pdbFileRef] = score

    #output = open("/media/bs674/ppi8t/NAMpopulation_alphafold/allNAMs/canonicalPorteins/deleteriousMutations/" + pan_group_id + "/usalignMatrix_TMScore", 'w')
    output = open("/home/songbs/deleteriousMutations/" + pan_group_id + "/usalignMatrix_TMScore", 'w')
    for row_label, row in zip(seq_names, distance_matrix):
        output.write('%s %s' % (row_label, ' '.join('%03s' % i for i in row)) + "\n")
    output.close()

    #output = open( "/media/bs674/ppi8t/NAMpopulation_alphafold/allNAMs/canonicalPorteins/deleteriousMutations/" + pan_group_id + "/usalignRMSD_Matrix", 'w')
    output = open( "/home/songbs/deleteriousMutations/" + pan_group_id + "/usalignRMSD_Matrix", 'w')
    for row_label, row in zip(seq_names, distance_matrix_RMSD):
        output.write('%s %s' % (row_label, ' '.join('%03s' % i for i in row)) + "\n")
    output.close()
