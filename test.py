def codonsToAmino(codons):
    aa_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 
        'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 
        'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 
        'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 
        'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 
        'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 
        'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 
        'TAA':'*', 'TAG':'*', 'TGA':'*',
        'TGC':'C', 'TGT':'C', 
        'TGG':'W',
    } #dictionary for Codon to AA nomenclature
    return aa_table.get(codons)

def aminoToCodon(amino):
    codon_table = {
    'A': ('GCT', 'GCC', 'GCA', 'GCG'),
    'C': ('TGT', 'TGC'),
    'D': ('GAT', 'GAC'),
    'E': ('GAA', 'GAG'),
    'F': ('TTT', 'TTC'),
    'G': ('GGT', 'GGC', 'GGA', 'GGG'),
    'I': ('ATT', 'ATC', 'ATA'),
    'H': ('CAT', 'CAC'),
    'K': ('AAA', 'AAG'),
    'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
    'M': ('ATG',),
    'N': ('AAT', 'AAC'),
    'P': ('CCT', 'CCC', 'CCA', 'CCG'),
    'Q': ('CAA', 'CAG'),
    'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
    'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
    'T': ('ACT', 'ACC', 'ACA', 'ACG'),
    'V': ('GTT', 'GTC', 'GTA', 'GTG'),
    'W': ('TGG',),
    'Y': ('TAT', 'TAC'),
    '*': ('TAA', 'TAG', 'TGA'),
    } #dictionary for AA to Codon nomenclature giving all options as tuples
    codon_table2 = {
    'A': 'GCT', 'A': 'GCC', 'A': 'GCA', 'A': 'GCG',
    'C': 'TGT', 'C': 'TGC',
    'D': 'GAT', 'D': 'GAC',
    'E': 'GAA', 'E': 'GAG',
    'F': 'TTT', 'F': 'TTC',
    'G': 'GGT', 'G': 'GGC', 'G': 'GGA', 'G': 'GGG',
    'I': 'ATT', 'I': 'ATC', 'I': 'ATA',
    'H': 'CAT', 'H': 'CAC',
    'K': 'AAA', 'K': 'AAG',
    'L': 'TTA', 'L': 'TTG', 'L': 'CTT', 'L': 'CTC', 'L': 'CTA', 'L': 'CTG',
    'M': 'ATG',
    'N': 'AAT', 'N': 'AAC',
    'P': 'CCT', 'P': 'CCC', 'P': 'CCA', 'P': 'CCG',
    'Q': 'CAA', 'Q': 'CAG',
    'R': 'CGT', 'R': 'CGC', 'R': 'CGA', 'R': 'CGG', 'R': 'AGA', 'R': 'AGG',
    'S': 'TCT', 'S': 'TCC', 'S': 'TCA', 'S': 'TCG', 'S': 'AGT', 'S': 'AGC',
    'T': 'ACT', 'T': 'ACC', 'T': 'ACA', 'T': 'ACG',
    'V': 'GTT', 'V': 'GTC', 'V': 'GTA', 'V': 'GTG',
    'W': 'TGG',
    'Y': 'TAT', 'Y': 'TAC',
    '*': 'TAA', '*': 'TAG', '*': 'TGA',
    } #dictionary for AA to Codon nomenclature giving single option as string
    return codon_table2.get(amino)

def DNAtoRNA(seq):
    table = {
    'A' : 'U', 
    'C' : 'G',
    'G' : 'C',
    'T' : 'A'
    }
    return table.get(seq)

def RNAtoDNA(seq):
    table = {
    'U' : 'A', 
    'G' : 'C',
    'C' : 'G',
    'A' : 'T'
    }
    return table.get(seq)

def info(): #main file for giving info about the specified fasta file
    file = open('read1.fasta', 'r') #file to read fasta from
    output = open('testoutput.out','w') #output file to give info on

    output.write("NOTE: The Amino Translations are only probable equivalent sequences due to the nature of Amino Acids having multiple equivalent DNA triplets, will try to fix later maybe\n\n")

    count = -1 #number of sequences in the file
    TypSeq = "NULL" #boolean for signalling whether passed sequence is an Amino Sequence or a DNA sequence

    for line in file: #each line is each "block" of genetic information, beginning with @NS, the for loop "reads" the first line
        TypSeq = "Unrecognized"
        count += 1 #counts how many sequences there are within the file

        seq = file.readline().strip('\n') #even though its first call of readline, reads the second line of each block due to for loop eating the first line
        if( (seq.count("A")+seq.count("G")+seq.count("C")+seq.count("T")) < len(seq) 
            and (seq.count("U")+seq.count("G")+seq.count("C")+seq.count("A")) < len(seq) 
            or  len(seq) <= 3 ): #determines type of sequence the input sequence is
            TypSeq = "Amino"
        elif("U" in seq):
            TypSeq = "RNA"
        else:
            TypSeq = "DNA"
        
        SeqDNA = "" #DNA sequence given the input sequence
        SeqAmino = "" #Amino sequence given the input sequence
        SeqRNA = "" #RNA sequence given the input sequence
        GC = 0 #GC % within each sequence
        if(TypSeq == "Amino"): 
            SeqAmino = seq
            
            for base in seq:
                SeqDNA = SeqDNA + str(aminoToCodon(base))

            for base in SeqDNA:
                SeqRNA = SeqRNA + str(DNAtoRNA(base))

            gCount = SeqDNA.count('G')
            cCount = SeqDNA.count('C')
            GC = (gCount + cCount) / (len(seq)*3) * 100

        elif(TypSeq == "DNA"):
            SeqDNA = seq
            
            for x in range( int(len(seq)/3) ):
                codons = str(seq[x*3]) + str(seq[x*3+1]) + str(seq[x*3+2])
                SeqAmino = SeqAmino + str(codonsToAmino(codons))
            
            for base in SeqDNA:
                SeqRNA = SeqRNA + str(DNAtoRNA(base))

            gCount = seq.count('G')
            cCount = seq.count('C')
            GC = (gCount + cCount) / (len(seq)) * 100

        elif(TypSeq == "RNA"):
            SeqRNA = seq

            for base in seq:
                SeqDNA = SeqDNA + str(RNAtoDNA(base))

            for x in range( int(len(SeqDNA)/3) ):
                codons = str(SeqDNA[x*3]) + str(SeqDNA[x*3+1]) + str(SeqDNA[x*3+2])
                SeqAmino = SeqAmino + str(codonsToAmino(codons))

            gCount = SeqDNA.count('G')
            cCount = SeqDNA.count('C')
            GC = (gCount + cCount) / (len(SeqDNA)) * 100

        output.write("[" + TypSeq + "] Sequence " + str(count) + ": " + str(seq) + "\n\tEquivalent DNA Sequence: " + str(SeqDNA) 
        + "\n\tEquivalent RNA Sequence: " + str(SeqRNA)
        + "\n\tEquivalent Amino Sequence: "+ str(SeqAmino) + "\n\tGC %: " + str(GC) + "\n\n") #write results to outfile

    file.close()
    output.close()

    print(count+1)

info()