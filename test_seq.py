from Bio.Seq import Seq

random_dna = Seq("attgcgcgtaaaatactgcccaccacgacggtctatgtcaccattgcgaa")
print(random_dna)
print(random_dna.complement())
print(random_dna.reverse_complement())