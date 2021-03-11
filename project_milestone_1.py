# question 1
def dna_count(dna):
    dna = dna.upper()
    count_A = dna.count('A')
    count_C = dna.count('C')
    count_G = dna.count('G')
    count_T = dna.count('T')
    return count_A, count_C, count_G, count_T

# question 2
def dna2rna(dna):
    rna = dna.replace('T', 'U')
    return rna.strip()

# question 3
def reverse_complement(dna):
    complement = ''
    for symbol in dna:
        if symbol == 'A':
            complement = complement + 'T'
        elif symbol == 'T':
            complement = complement + 'A'
        elif symbol == 'C':
            complement = complement + 'G'
        else:
            complement = complement +'C'
    reverse = complement[::-1] #this reverses the string
    return reverse
dna = "AAAACCCGGT"
print("reversed complement is:", reverse_complement(dna))

# question 4
def mendels_law (hom, het, rec):
    pop_total = hom + het + rec
    ress = (het/pop_total) * (rec - 1)/(pop_total - 1)
    het_ress = (het/pop_total) * het/(pop_total - 1) * 0.5
    ress_het = (rec/pop_total) * het/(pop_total - 1) * 0.5
    het_het = het/pop_total * (het - 1)/(pop_total -1) * 0.25

    return (1 - (ress + het_ress + ress_het + het_het))

print (mendels_law (2, 2, 2))

# question 5
def fibonacci_rabbits(n ,k):
    if n == 0:
        return 0
    if n == 1:
        return 1
    else:
        return fibonacci_rabbits(n - 1, k) + k * fibonacci_rabbits(n - 2, k)
print (fibonacci_rabbits (5, 3))

# question 6
def GC_content(dna_list):
    count_GC = []
    for i in range(len(dna_list)):
        dna_str = dna_list[i]
        count_GC.append(dna_str.count('G') + dna_str.count('C'))
    maxGC_num = max(count_GC)
    maxGC_index = count_GC.index(maxGC_num)
    perc_GC = ((maxGC_num)/(len(dna_list(maxGC_index))))*100
    return (maxGC_index, perc_GC)

# question 7
def rna2codon(seq):
    table = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',

        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',

        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    protein = ""
    if len(seq)%3 ==0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i +3]
            if table[codon]!="STOP":
                    protein += table[codon]
    return protein
print( rna2codon("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"))

# question 8
def locate_substring(dna_snippet, dna):
        indexes = [i for i in range(len(dna_snippet)) if dna_snippet.startswith(dna, i)]
        return indexes
print(locate_substring("GATATATGCATATACTT","ATAT"))

# question 9
def hamming_dist(dna1, dna2):
    count = 0
    i = 0
    while( i < len(dna1)):
        if dna1[i] != dna2[i]:
            count += 1
        i += 1
    return count
a = 'GAGCCTACTAACGGGAT'
b = 'CATCGTAATGACGGCCT'
print(hamming_dist(a,b))

# question 10
def count_dom_phenotype(genotypes):
    AAAA=genotypes[0]*2
    AAAa=genotypes[1]*2
    AAaa=genotypes[2]*2
    AaAa=genotypes[3]*1.5
    Aaaa=genotypes[4]*0.5
    aaaa=genotypes[5]*0
    answer= AAAA+AAAa+AAaa+AaAa+Aaaa+aaaa
    return answer

# question 11
def source_rna(protein):
    prosteam=protein +'*'
    if '*' in protein:
        prosteam=protein
    triplebase= {
        'F':2,
        'L':6,
        'I':3,
        'M':1,
        'V':4,
        'S':6,
        'P':4,
        'T':4,
        'A':4,
        'Y':2,
        'H':2,
        'Q':2,
        'N':2,
        'K':2,
        'D':2,
        'E':2,
        'C':2,
        'W':1,
        'R':6,
        'G':4,
        '*':3   }
    total=1
    for wow in prosteam:
        total*=triplebase[wow]
    if len(protein.strip())>0:
        return total%1000000
    else:
        return 0

# question 12
def splice_rna(dna, intron_list):
    yra = dna
    for intron in intron_list:
        yra = yra.replace(intron,"")
    return rna2codons(dna2rna(yra))
