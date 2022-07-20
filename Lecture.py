
def revSeq(seq):
    data = seq[::-1]
    compDic = {"A":"T","T":"A","G":"C","C":"G"}
    result = ""
    for i in data:
        result += compDic[i]
    return result

def DNAtoRNA(revComp_result):
    RNA_seq = revComp_result.replace("T","U")
    return RNA_seq

def RNAtoProtein(RNA_result):
    protein = ""
    Codon = {'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C', 'UUC': 'F', 'UCC': 'S', 
    'UAC': 'Y', 'UGC': 'C','UUA': 'L', 'UCA': 'S', 'UAA': 'STOP', 'UGA': 'STOP','UUG': 'L', 
    'UCG': 'S', 'UAG': 'STOP', 'UGG': 'W','CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R', 
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R','CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 
    'CGA': 'R', 'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R','AUU': 'I', 'ACU': 'U', 
    'AAU': 'N', 'AGU': 'S', 'AUC': 'I', 'ACC': 'U', 'AAC': 'N', 'AGC': 'S','AUA': 'I', 
    'ACA': 'U', 'AAA': 'K', 'AGA': 'R', 'AUG': 'M', 'ACG': 'U', 'AAG': 'K', 'AGG': 'R',
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G', 'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 
    'GGC': 'G','GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'} 
    end = len(RNA_result) - (len(RNA_result)%3)
    for i in range(0,end,3):
        codon = RNA_result[i:i+3]
        if Codon[codon] == 'STOP': protein += "None";break
        protein += Codon[codon]
    return protein 


def main():
    seq = input("Enter the sequence information!!! : ")
    print("################################################")
    print("### Input Sequence information :",seq)
    print("################################################")
    revComp_result = revSeq(seq)
    print("################################################")
    print("### Complementary Sequence information :",revComp_result)
    print("################################################")
    RNA_result = DNAtoRNA(revComp_result)
    print("################################################")
    print("### RNA Sequence information :",RNA_result)
    print("################################################")
    Protein_result = RNAtoProtein(RNA_result)
    print("################################################")
    print("### RNA to Protein information :",Protein_result)
    print("################################################")

    File_name = input("Plz Enter the save file name here : ")
    File_name2 = str(File_name) + ".txt"
    file_result = open(File_name2,"w")
    file_result.write("################################################"+ "\n")
    file_result.write("Input sequence information : " + seq + "\n")
    file_result.write("################################################"+ "\n")
    file_result.write("Complementary Sequence information : " + revComp_result + "\n")
    file_result.write("################################################"+"\n")
    file_result.write("RNA Sequence information : " + RNA_result + "\n")
    file_result.write("################################################"+ "\n")
    file_result.write("Protein Sequence information : " + Protein_result + "\n")
    file_result.write("################################################"+ "\n")
    file_result.close()
    print("Result saved on file to : ", File_name2)
    

main()