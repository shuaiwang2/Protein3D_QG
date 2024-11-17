import os
import glob

def read_mafft_alignment(file_path):
    alignment = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                sequence_name = line.strip()[1:]
                alignment[sequence_name] = ''
            else:
                alignment[sequence_name] += line.strip()
    return alignment

def count_differences(seq1, seq2):
    differences1 = 0
    differences2 = 0
    for base1, base2 in zip(seq1, seq2):
        if base1 != base2:
            if base1 == '-' or base2 == '-':
                differences1 += 1
            else:
                differences2 += 1
    return differences1,differences2

    
def main():
    directory = '/media/shuaiwang/3d/NAMpopulation_alphafold/allNAMs/canonicalPorteins/deleteriousMutations'
    folders = [folder for folder in os.listdir(directory) if os.path.isdir(os.path.join(directory, folder))]
    #folders = [folder for folder in os.listdir(directory) if os.path.isdir(os.path.join(directory, folder)) and folder.startswith('200')]

    data_dict = {}
    for folder_name in folders:
        file_paths = glob.glob(os.path.join(directory, folder_name, 'sequences.mafft'))
        for file_path in file_paths:
            #with open(file_path) as mafft:
            #    for pdbline in mafft:
                #    if (pdbline.startswith('ATOM')) and pdbline[13:15] == 'CA':
                #        print("B73",folder_name,int(pdbline[22:26].strip()),float(pdbline[60:66]),sep='\t') 
                #mafft_file = "sequences.mafft"
                alignment = read_mafft_alignment(file_path)
                result = 0
                difference_res =0
                sequence_names = list(alignment.keys())
                for i in range(len(sequence_names)):
                    for j in range(i+1, len(sequence_names)):
                        sequence1 = alignment[sequence_names[i]]
                        sequence2 = alignment[sequence_names[j]]
                        differences1,differences2 = count_differences(sequence1, sequence2)
                        #print(folder_name,sequence_names[i],sequence_names[j],differences,differences/len(alignment[sequence_names[0]]),sep='\t')
                        print(folder_name,sequence_names[i],sequence_names[j],differences1,differences2,sep ='\t')
                        #difference_res+= differences
                        #if differences == 0:
                        #    result +=1
                #print(folder_name,difference_res/(len(sequence_names)- 1)/len(sequence_names)/len(alignment[sequence_names[0]]),result/(len(sequence_names)- 1)/len(sequence_names)*2,len(sequence_names),sep='\t')
if __name__ == "__main__":
    main()
