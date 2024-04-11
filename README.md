# NLR_Sordariales
Script used in the incoming article "Systematic exploration of domain assortments in NOD-like receptors uncovers two types of NACHT domains in Sordariales fungi" by Bonometti L. et al.

1) script_correction_exonerate_v10.py
This script was designed to analyze EXONERATE outputs from the alignment of NLR sequences against a single whole genome and the alignment of their NB domain sequences against the same genome. It also analyzes the original gene prediction. Based on the alignment of homologous genes on mispredicted or unpredicted genes, it :
  - corrects genes that were erroneously predicted as two or more different genes;
  - prolongates genes that were erroneously truncated in C-term or N-term;
  - prolongates genes that do not start with methionine or do not end with stop codon;
  - predicts new genes that were missed by original gene prediction;
  - corrects gene sequences when original gene prediction includes stop codons inside the sequence;
  - identifies pseudo-genes when an NLR gene cannot be predicted without stop codons in its sequence.
The script takes as input:
  - the fasta containing all the NLRs, as nucleotide sequences (optional);
  - a directory containing the nucleotide fasta files of all the NLR orthogroups, with both the NLRs and the non-NLRs genes of those orthogroups (optional);
  - the GFF file of the original gene prediction associated with the assembled genome (optional);
  - the assembly fasta file of the genome (mandatory);
  - the file containing the exonerate outputs for the alignment of all NB sequences on the analyzed genome, without gff included (mandatory);
  - the threshold for accepting an exonerate hit of an NLR as a valid hit (optional, default = 100);
  - the threshold for accepting an exonerate hit of an NB as a valid hit (optional, default = 100);
  - the minimal size, in bp, for accepting a valid NLR (optional, default = 1000);
  - the proportion, in percent, of a query gene that can be shortened in the C-terminal extremity of the target hit (optional, default = 10);
The script outputs:
  - a fasta file of all the NLRs finally identified in the genome, in nucleotides;
  - a fasta file of all the NLRs finally identified in the genome, in amino acids;
  - a text file listing the identified pseudo-genes;
  - a text file providing the list of NLR genes whose prediction did not change from the original gene prediction, the list of NLR genes whose prediction was modified from the original gene prediction, with their exons, the list of NLR genes not yet identified as NLRs, with their exons, the list of NLR gene that were missed by the original gene prediction, and the list of identified pseudo-genes.
  - a text file summarizing the analysis and its results.
