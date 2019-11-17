# protein-align

Script for creating matrices of substitutions, deletions, and insertions of amino acids from SABmark sup alignments.

gen_nseq() 
PURPOSE: Finds the number of distinct sequences in all groups
INPUT: none
OUTPUT: dictionary with {group number : nseq}
Path directory dependent

neighbor(seq, i)
PURPOSE: Finding neighbors of a protein at index i in sequence seq; identifying sub/in/del
INPUT: seq - string ; i - integer
OUTPUT: list [left neighbor, current protein, right neighbor, length]
For substitutions, length = 1
For indels, length = length of indel, and current protein = "-"
For proteins neighboring indels, length = length of indel + 1, current protein != "-"

get_seqs(file)
PURPOSE: Extract sequences from a .fasta file
INPUT: .fasta file
OUTPUT: array of two strings (sequences), with seq 0 being the ancestor/parent, and seq 1 being the descendant/child
in_matrix(seqs) 	#INCOMPLETE
PURPOSE: Track frequencies of insertions in one alignment (2 seq)
INPUT: array of two string sequences, with the 0th being the parent and the 1st the child
OUTPUT: Dictionary {key : count} ; key is a string of len 3-4 in the format [parent leftFlank, parent rightFlank, length of insert]
dictionary[‘type’] will get you ‘in’
Use in conjunction with get_seqs

del_matrix
PURPOSE: Track frequencies of deletions in one alignment (2 seq)
INPUT: array of two string sequences, with the 0th being the parent and the 1st the child
OUTPUT: Dictionary {key : count} ; key is a string of len 3-4 in the format [parent leftFlank, parent rightFlank, length of deletion]
dictionary[‘type’] will get you ‘del’
sub_matrix 		#NEEDSUPDATING
PURPOSE: Track frequencies of insertions in one alignment (2 seq)
INPUT: array of two string sequences, with the 0th being the parent and the 1st the child
OUTPUT: Dictionary {key : count} ; key is a string of len 4 in the format [parent leftFlank, current protein, parent rightFlank, subbed protein]
dictionary[‘type’] will get you ‘sub’
