"""
DNA Sequence Analyzer

Combining transcription, translation, mutation, and ORF analysis.
"""

# Universal Genetic Code Table (RNA and DNA versions)
CODON_TABLE_RNA = {
    'AUA':'Ile', 'AUC':'Ile', 'AUU':'Ile', 'AUG':'Met',
    'ACA':'Thr', 'ACC':'Thr', 'ACG':'Thr', 'ACU':'Thr',
    'AAC':'Asn', 'AAU':'Asn', 'AAA':'Lys', 'AAG':'Lys',
    'AGC':'Ser', 'AGU':'Ser', 'AGA':'Arg', 'AGG':'Arg',                 
    'CUA':'Leu', 'CUC':'Leu', 'CUG':'Leu', 'CUU':'Leu',
    'CCA':'Pro', 'CCC':'Pro', 'CCG':'Pro', 'CCU':'Pro',
    'CAC':'His', 'CAU':'His', 'CAA':'Gln', 'CAG':'Gln',
    'CGA':'Arg', 'CGC':'Arg', 'CGG':'Arg', 'CGU':'Arg',
    'GUA':'Val', 'GUC':'Val', 'GUG':'Val', 'GUU':'Val',
    'GCA':'Ala', 'GCC':'Ala', 'GCG':'Ala', 'GCU':'Ala',
    'GAC':'Asp', 'GAU':'Asp', 'GAA':'Glu', 'GAG':'Glu',
    'GGA':'Gly', 'GGC':'Gly', 'GGG':'Gly', 'GGU':'Gly',
    'UCA':'Ser', 'UCC':'Ser', 'UCG':'Ser', 'UCU':'Ser',
    'UUC':'Phe', 'UUU':'Phe', 'UUA':'Leu', 'UUG':'Leu',
    'UAC':'Tyr', 'UAU':'Tyr', 'UAA':'Stop', 'UAG':'Stop',
    'UGC':'Cys', 'UGU':'Cys', 'UGA':'Stop', 'UGG':'Trp',
}

CODON_TABLE_DNA = {
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

def transcription(dna_template):
    """
    Converts a 3'-5' DNA template strand into its complementary 5'-3' mRNA sequence.
    Universal codes included for broader compatibility.
    """
    comp_map = str.maketrans('ACGTRYSWKMBDHVN', 'UGCAYRSWMKVHDBN')
    # Remove spaces/formatting
    dna_clean = dna_template.upper().replace(" ", "").replace("-", "")
    # Base pairing rules: A->U, T->A, C->G, G->C
    mrna = dna_clean.translate(comp_map)
    return mrna

def find_codons(sequence):
    """
    Divides a sequence into codons (groups of three nucleotides).
    """
    sequence = sequence.upper().replace(" ", "").replace("-", "")
    return [sequence[i:i+3] for i in range(0, len(sequence) - len(sequence) % 3, 3)]

def translate_codon(codon):
    """
    Translates a single codon (DNA or RNA) into an amino acid.
    """
    codon = codon.upper()
    if 'U' in codon:
        return CODON_TABLE_RNA.get(codon, '?')
    else:
        return CODON_TABLE_DNA.get(codon, '?')

def translate_sequence(sequence):
    """
    Translates a full mRNA or DNA sequence into an amino-acid string.
    """
    codons = find_codons(sequence)
    protein = []
    for codon in codons:
        aa = translate_codon(codon)
        protein.append(aa)
    return "-".join(protein)

def find_orf(dna_sequence):
    """
    Identifies a possible Open Reading Frame in a 5'-3' DNA coding sequence.
    """
    dna_sequence = dna_sequence.upper().replace(" ", "").replace("-", "")
    
    # Find start codon
    start_idx = dna_sequence.find("ATG")
    
    if start_idx == -1:
        return {
            "start_pos": None,
            "stop_pos": None,
            "codons": [],
            "orf_seq": "",
            "protein": "",
            "is_valid": False
        }
        
    # Extract candidate sequence starting from ATG
    candidate_seq = dna_sequence[start_idx:]
    
    # Divide into codons
    codons = find_codons(candidate_seq)
    
    # Find in-frame stop codon and translate
    stop_codons = {"TAA", "TAG", "TGA"}
    orf_codons = []
    protein = []
    stop_pos = None
    stop_codon_found = None
    
    for i, codon in enumerate(codons):
        orf_codons.append(codon)
        if codon in stop_codons:
            stop_pos = start_idx + (i * 3)
            stop_codon_found = codon
            break
        aa = translate_codon(codon)
        protein.append(aa)
        
    is_valid = stop_pos is not None
    
    return {
        "start_pos": start_idx,
        "stop_pos": stop_pos,
        "stop_codon": stop_codon_found,
        "codons": orf_codons,
        "orf_seq": "".join(orf_codons),
        "protein": "-".join(protein),
        "is_valid": is_valid
    }

if __name__ == "__main__":
    print("=== DNA Sequence Analyzer ===")
    
    # Accept user input as requested in final challenge
    try:
        dna = input("Enter a DNA sequence (defaulting to Activity 2 seq if empty): ").strip()
    except EOFError:
        dna = ""
        
    if not dna:
        dna = "CCGATGAAACCGTACGGGTAACTG"
        
    print(f"\nDNA sequence: {dna}")
    
    mrna = dna.upper().replace("T", "U")
    print(f"mRNA sequence: {mrna}")
    
    codons = find_codons(dna)
    print(f"Codons: {' | '.join(codons)}")
    
    protein = translate_sequence(dna)
    print(f"Protein sequence: {protein}")
    
    orf_data = find_orf(dna)
    print(f"Start codon position: {orf_data['start_pos']}")
    print(f"Stop codon position: {orf_data['stop_pos']}")
    print(f"ORF sequence: {orf_data['orf_seq']}")
    print(f"ORF status (valid): {orf_data['is_valid']}")

    print("\n=== Activity 1: Decoding and Mutation Challenge ===")
    orig_template = "TAC CGA TTT ACC ACT"
    mut_template = "TAC CGA TAT ACC ACT"
    
    orig_mrna = transcription(orig_template)
    mut_mrna = transcription(mut_template)
    
    orig_prot = translate_sequence(orig_mrna)
    mut_prot = translate_sequence(mut_mrna)
    
    print(f"Original template: {orig_template}")
    print(f"Original mRNA: {orig_mrna}")
    print(f"Original protein: {orig_prot}")
    
    print(f"Mutated template: {mut_template}")
    print(f"Mutated mRNA: {mut_mrna}")
    print(f"Mutated protein: {mut_prot}")
    
    print(f"Mutation changes protein: {orig_prot != mut_prot}")

    print("\n=== Activity 2: ORF Challenge ===")
    orig_orf = "CCGATGAAACCGTACGGGTAACTG"
    mut_orf = "CCGATGAAACCGTACGGGTGAACTG"
    
    res1 = find_orf(orig_orf)
    res2 = find_orf(mut_orf)
    
    print(f"Original sequence: {orig_orf}")
    print(f"Start: {res1['start_pos']}")
    print(f"Stop: {res1['stop_pos']} ({res1['stop_codon']})")
    print(f"Protein: {res1['protein']}")
    
    print(f"\nModified sequence: {mut_orf}")
    print(f"Start: {res2['start_pos']}")
    print(f"Stop: {res2['stop_pos']} ({res2['stop_codon']})")
    print(f"Protein: {res2['protein']}")
    
    print(f"\nEffect of mutation: Stop codon is lost/shifted, making ORF invalid" if not res2['is_valid'] else f"\nEffect of mutation: Protein changed to {res2['protein']}")
