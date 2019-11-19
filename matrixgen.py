import glob
from pathlib import Path
import csv
from aminoacids import *

###############
# DESCRIPTION #
###############
"""
This file serves to generate three matrices based off of the SABmark
sequence alignments. These three matrices are the substitution, insertion,
and deletion matrices, which can be represented in dictionary format
or written to a CSV file.

Author Yolanda Shen
"""


# Variables to be modified

matrix_type = 'in'

generate_matrix = False
print_matrix_stats = False
generate_csv = False

first_group = 1
last_group = 425

path_to_sup = 'C:/Users/yolan/SABmark/SABmark/sup/'
path_to_csv_file = path_to_sup + 'in_matrix.csv'


#############################################
# Do not modify anything below this comment #
#############################################

matrix = {}
files_parsed = 0
expected_parsed = 0

# Generates a dictionary containing then number of distinct sequences in each group
# {group number : nseq}
def gen_nseq():
    nseq = {}
    setSummary = Path(path_to_sup).glob('set.summary').__iter__()
    summaryIter = open(setSummary.__next__())
    summaryIter.__next__()
    for line in summaryIter:
        line = line.split('\t')[:2]
        group = line[0]
        num = line[1]
        if group not in nseq.keys():
            nseq[int(group)] = int(num)
    return nseq

# Returns neighbors of the protein at index i
# If i is an indel ("-"), returns the nearest non-indel neighbor
# OUTPUT [left neighbor, current protein, right neighbor, length]
# '$' indicates beginning/end of sequence
# For substitutions, length = 1
# For indels, length = length of indel, and current protein = "-"
# For proteins neighboring indels, length = length of indel + 1, 
#   current protein != "-"
def neighbor(seq, i):
    count = 1
    def left(seq, i):
        nonlocal count
        if i == 0:
            return "$"
        elif seq[i-1] == "-":
            count += 1
            return left(seq, i-1)
        else:
            return seq[i-1]
    def right(seq, i):
        nonlocal count
        if i == len(seq)-1:
            return "$"
        elif seq[i+1] == "-":
            count += 1
            return right(seq, i+1)
        else:
            return seq[i+1]
    return [left(seq, i), seq[i], right(seq, i), count]

# Gets the protein sequences from a file, 
# returns as an array of two strings.
def get_seqs(file):
    iter = open(file).__iter__()
    line = iter.__next__()
    seqs = []
    seq = ""
    while iter:
        if ">" in line:
            if seq != "":
                seqs.append(seq)
                seq = ""
            line = iter.__next__()
        else:
            seq += line[:len(line)-1]
            try:
                line = iter.__next__()
            except:
                StopIteration
                seqs.append(seq)
                break
    if len(seq[0]) != len(seq[1]):
        return Exception("Alignment sequences different lengths")
    return seqs

# Insertion matrix on a given fasta alignment 
# Returns a dictionary of frequencies with 3-4 char key
# letters 0-1 = parent l, r
# chars 2-3 = number of inserted proteins

# insertion is where p_char == '-' and c_char != '-'
# insertion length is dependent on parent[3]
# upon finding an insertion, make sure that it doesn't duplicate count
def in_matrix(seqs):
    parent = seqs[0]
    child = seqs[1]
    assert(len(parent) == len(child))
    global matrix
    repeat = False
    for i in range(len(parent)):
        p_neighbor = neighbor(parent, i)
        c_neighbor = neighbor(child, i)
        if p_neighbor[1] == "-" and not repeat:
            repeat = True
            j = i
            while j < len(parent)-1 and parent[j+1] == '-':
                j += 1
            seq = p_neighbor[0] + neighbor(parent, j)[2] + str(p_neighbor[3])
            if seq in matrix.keys():
                matrix[seq] += 1
            else:
                matrix[seq] = 1
        else:
            repeat = parent[i] == "-"
    return matrix

# Deletion matrix on a given fasta alignment
# Returns a dictionary of frequencies with 3-4 char key
# letters 0-1 = parent l, r
# chars 2-3 = number of deleted proteins

# deletion is where p_char != '-' and c_char == '-'
def del_matrix(seqs):
    parent = seqs[0]
    child = seqs[1]
    assert(len(parent) == len(child))
    global matrix
    repeat = False
    for i in range(len(parent)):
        p_neighbor = neighbor(parent, i)
        c_neighbor = neighbor(child, i)
        if c_neighbor[1] == "-" and not repeat:
            repeat = True
            j = i
            while j < len(parent)-1 and child[j+1] == '-':
                j += 1
            seq = p_neighbor[0] + neighbor(parent, j)[2] + str(c_neighbor[3])
            if seq in matrix.keys():
                matrix[seq] += 1
            else:
                matrix[seq] = 1
        else:
            repeat = child[i] == "-"
    return matrix

# Substitution matrix on a given fasta alignment (2 seq)
# Returns a dictionary of frequencies with 4 letter key
# letters 0-2 = parent left, curr, right
# letter 3 = child sub at sam position as curr
def sub_matrix(seqs):
    parent = seqs[0]
    child = seqs[1]
    assert(len(parent) == len(child))
    global matrix
    for i in range(len(parent)):
        p_neighbor = neighbor(parent, i)
        c_neighbor = neighbor(child, i)
        if p_neighbor[1] != c_neighbor[1] \
            and p_neighbor[3] == c_neighbor[3] == 1 \
                and c_neighbor[1] != '-' \
                    and p_neighbor[1] != "-":
            seq = "".join(p_neighbor[:3] + c_neighbor[1:2])
            if seq in matrix.keys():
                matrix[seq] += 1
            else:
                matrix[seq] = 1
    return matrix

# Calls the corresponding function on a file, given the matrix type
def apply_matrix(str_type, file):
    global matrix
    funcs = {"in": in_matrix, "del": del_matrix, "sub": sub_matrix}
    funcs[str_type](get_seqs(file))
    return

# Iterates through the files in a group (does not count repeats, 
# and the parent sequence is consistent)
def file_iterator(str_type, path, group):
    global matrix
    count = 0
    gen = path.__iter__()
    while gen and count < nseq[group]:
        try:
            apply_matrix(str_type, gen.__next__())
            count += 1
        except:
            StopIteration
            break
    return count

# Iterates through specified group range to generate the matrix
def matrix_generator(str_type, bool):
    if not bool:
        return
    global matrix
    global files_parsed
    global expected_parsed
    if len(matrix) == 0 or matrix["type"] != str_type:
        matrix = {}
        matrix["type"] = str_type
    for group in range(first_group, last_group + 1):
        print("Analyzing... %d/%d complete"%(group, last_group))
        files_parsed += file_iterator(str_type, glob.iglob(path_to_sup + 'group' + str(group) + '/reference/*.fasta'), group)
        expected_parsed += nseq[group]

# Writes the generated matrix to a CSV file
def write_to_csv(csv_file, bool):
    if not bool:
        return
    with open(csv_file, 'w') as f:
        for key in matrix.keys():
            if key == 'type':
                if matrix[key] == 'sub':
                    f.write("Lflank,Mid,Rflank,Sub,Count")
                elif matrix[key] == 'del':
                    f.write("Lflank,Rflank,Del Len,Count")
                elif matrix[key] == 'in':
                    f.write("Lflank,Rflank,In Len,Count")
            elif matrix['type'] == 'sub':
                f.write("%s,%s,%s,%s,%s\n" % (key[0], key[1], key[2], key[3], matrix[key]))
            else:
                f.write("%s,%s,%s,%s\n" % (key[0], key[1], key[2:], matrix[key]))
    print("CSV successfully written.")

# Prints the type, length, and count of fasta alignments analyzed of the matrix
def print_matrix_info(bool):
    if not bool:
        return
    print("Matrix type: " + matrix['type'])
    print("Length of matrix: " + str(len(matrix) - 1))
    print("Number of fasta alignments parsed: " + str(files_parsed))
    print("Expected number of files parsed: " + str(expected_parsed))

# These do the function calling 
nseq = gen_nseq()
matrix_generator(matrix_type, generate_matrix)
print_matrix_info(print_matrix_stats)
write_to_csv(path_to_csv_file, generate_csv)
