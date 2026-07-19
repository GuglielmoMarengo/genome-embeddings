class Genome:
    def __init__(self, sequence, header=None):
        self.validate_string("Sequence", sequence)

        sequence = sequence.upper()
        valid_nucleotides = {"A", "C", "G", "T"}

        for i, character in enumerate(sequence, start=1):
            if character not in valid_nucleotides:
                raise ValueError(f"Invalid character {character} in sequence at position {i}.")

        self.sequence = sequence
        self.header = header

    @classmethod
    def from_fasta(cls, filepath):
        with open(filepath) as file:
            lines = file.readlines()

            if not lines:
                raise ValueError("FASTA file is empty.")

            header = lines[0].strip()
            if not header.startswith(">"):
                raise ValueError("FASTA header must start with '>'.")

            sequence_lines = lines[1:]
            sequence = "".join(line.strip() for line in sequence_lines)

            return cls(sequence=sequence, header=header)

    @staticmethod
    def validate_string(name, value):
        if type(value) != str:
            raise TypeError(f"{name} must be a string.")
        if len(value) == 0:
            raise ValueError(f"{name} cannot be empty.")

    def length(self):
        return len(self.sequence)

    def gc_content(self):
        gc_bases = {"C", "G"}
        gc_count = 0

        for character in self.sequence:
            if character in gc_bases:
                gc_count += 1
        
        return gc_count / self.length()

    def reverse_complement(self):
        reverse_complement = []

        reversed_sequence = self.sequence[::-1]

        complement = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
        } 
        
        for character in reversed_sequence:
            reverse_complement.append(complement[character])
        
        return "".join(reverse_complement)
    
    def kmers(self, k):
        if type(k) != int:
            raise TypeError(f"k must be an integer, got {type(k).__name__}.")
        if k <= 0:
            raise ValueError(f"{k} must be positive.")
        if k > self.length():
            raise ValueError(f"{k} k cannot exceed the sequence length.")
        
        kmers = []
        for i in range(self.length() - k + 1):
            kmer = self.sequence[i:i+k]
            kmers.append(kmer)
        return kmers
    
    def kmer_frequencies(self, k):
        if type(k) != int:
            raise TypeError(f"k must be an integer, got {type(k).__name__}.")
        if k <= 0:
            raise ValueError(f"{k} must be positive.")
        if k > self.length():
            raise ValueError(f"{k} cannot exceed the sequence length.")

        frequencies = {}

        for kmer in self.kmers(k):
            if kmer in frequencies:
                frequencies[kmer] += 1
            else:
                frequencies[kmer] = 1
        
        return frequencies