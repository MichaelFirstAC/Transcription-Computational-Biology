"""
DNA Sequence Toolkit
Combines validation, frequency analysis, reverse complement, and translation.
"""

# 1. Global Constants
DNA_CODON_TABLE = {
    'ATA':'Ile', 'ATC':'Ile', 'ATT':'Ile', 'ATG':'Met',
    'ACA':'Thr', 'ACC':'Thr', 'ACG':'Thr', 'ACT':'Thr',
    'AAC':'Asn', 'AAT':'Asn', 'AAA':'Lys', 'AAG':'Lys',
    'AGC':'Ser', 'AGT':'Ser', 'AGA':'Arg', 'AGG':'Arg',                 
    'CTA':'Leu', 'CTC':'Leu', 'CTG':'Leu', 'CTT':'Leu',
    'CCA':'Pro', 'CCC':'Pro', 'CCG':'Pro', 'CCT':'Pro',
    'CAC':'His', 'CAT':'His', 'CAA':'Gln', 'CAG':'Gln',
    'CGA':'Arg', 'CGC':'Arg', 'CGG':'Arg', 'CGT':'Arg',
    'GTA':'Val', 'GTC':'Val', 'GTG':'Val', 'GTT':'Val',
    'GCA':'Ala', 'GCC':'Ala', 'GCG':'Ala', 'GCT':'Ala',
    'GAC':'Asp', 'GAT':'Asp', 'GAA':'Glu', 'GAG':'Glu',
    'GGA':'Gly', 'GGC':'Gly', 'GGG':'Gly', 'GGT':'Gly',
    'TCA':'Ser', 'TCC':'Ser', 'TCG':'Ser', 'TCT':'Ser',
    'TTC':'Phe', 'TTT':'Phe', 'TTA':'Leu', 'TTG':'Leu',
    'TAC':'Tyr', 'TAT':'Tyr', 'TAA':'Stop', 'TAG':'Stop',
    'TGC':'Cys', 'TGT':'Cys', 'TGA':'Stop', 'TGG':'Trp',
}

# 2. Validation
def validate_dna(dna_seq):
    """ Checks if DNA sequence contains only valid IUPAC nucleotide codes. """
    valid_iupac = set("ACGTURYKMSWBDHVN")
    return set(dna_seq.upper()).issubset(valid_iupac)

# 3. Frequency & Ambiguity Analysis
def frequency_advanced(seq):
    """ Calculates symbol frequency, ACGT percentages, and checks for ambiguous codes. """
    dic = {}
    seq_upper = seq.upper()
    total_len = len(seq_upper)
    
    for s in seq_upper:
        dic[s] = dic.get(s, 0) + 1
        
    percentages = {base: (dic.get(base, 0) / total_len * 100) if total_len > 0 else 0 
                   for base in ['A', 'C', 'G', 'T']}
    
    standard_bases = {'A', 'C', 'G', 'T'}
    has_ambiguous = not set(dic.keys()).issubset(standard_bases)
    
    return {
        "frequency": dic,
        "acgt_percentages": percentages,
        "contains_ambiguous": has_ambiguous
    }

# 4. Reverse Complement
def reverse_complement_optimized(dna_seq):
    """ Computes the reverse complement of the DNA sequence efficiently, including all IUPAC codes. """
    if not validate_dna(dna_seq):
        raise ValueError("Invalid DNA sequence")
        
    comp_map = str.maketrans('ACGTRYSWKMBDHVNU', 'TGCAYRSWMKVHDBNA')
    return dna_seq.upper().translate(comp_map)[::-1]

# 5. Translation
def translate_codon(cod):
    """ Translates a single DNA codon into an amino acid. """
    return DNA_CODON_TABLE.get(cod.upper())

def translate_sequence(dna_seq):
    """ Translates a full DNA sequence into an amino acid string. """
    if not validate_dna(dna_seq):
        raise ValueError("Invalid DNA sequence")
    
    dna_seq = dna_seq.upper()
    protein = []
    
    # Loop through the sequence in chunks of 3 (ignoring trailing bases that don't make a full codon)
    for i in range(0, len(dna_seq) - len(dna_seq) % 3, 3):
        codon = dna_seq[i:i+3]
        amino_acid = translate_codon(codon)
        
        if amino_acid:
            protein.append(amino_acid)
            if amino_acid == 'Stop':  # Stop codon encountered
                break
                
    return "-".join(protein)

# --- Test Execution Block ---
if __name__ == "__main__":
    # A test sequence that starts with a start codon (ATG) and ends with a stop codon (TAA)
    test_seq = "ATGCGATACGCTTACGGGTAA"
    
    print(f"Original Sequence: {test_seq}")
    print(f"Is Valid DNA?: {validate_dna(test_seq)}")
    
    stats = frequency_advanced(test_seq)
    print(f"Contains Ambiguous Codes: {stats['contains_ambiguous']}")
    print(f"Base Frequencies: {stats['frequency']}")
    
    rev_comp = reverse_complement_optimized(test_seq)
    print(f"Reverse Complement: {rev_comp}")
    
    protein = translate_sequence(test_seq)
    print(f"Translated Protein: {protein}")