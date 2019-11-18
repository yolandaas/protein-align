from math import *
import csv
from aminoacids import *

###############
# DESCRIPTION #
###############
"""
This file serves to generate a substitution frequency matrix
based on the BLOSUM62 model. The equations used are identical to
those found on the BLOSUM62 NCBI file and BLOSUM Wikipedia page.
The file requires import from a csv file that includes these columns:

Mid - amino acid i, which is substituted by j
Sub - amino acid j,which is substituted for i
Count - count of observed substitutions, given i, j

Author Yolanda Shen

IN PROGRESS
"""

# Variables to be modified
sub_csv = 'C:/Users/yolan/SABmark/SABmark/sup/sub_matrix.csv'
csv_file = 'C:/Users/yolan/SABmark/SABmark/sup/subfreq_matrix.csv'
generate_csv = True

#############################################
# Do not modify anything below this comment #
#############################################

# Returns the log odd ratio of a substitution
def log_odd_ratio(obs_prob, exp_prob):
    return 2 * log2(obs_prob / exp_prob)

# Function to calculate the matrix values for BLOSUM62
def score_equation(scale_factor, prob_ij, qi, qj):
    return (1 / scale_factor) * log(prob_ij / (qi * qj))

# Finds the column numbers corresponding to aforementioned
# column headers in the csv file
def column_headers():
    global mid_col
    global sub_col
    global count_col
    global line_count
    assert line_count == 0
    headers = csv_reader.__next__()
    line_count += 1
    for i in range(len(headers)):
        header = headers[i]
        if header == 'Mid':
            assert mid_col == -1, "multiple Mid columns"
            mid_col = i
        elif header == 'Sub':
            assert sub_col == -1, "multiple Sub columns"
            sub_col = i
        elif header == 'Count':
            assert count_col == -1, "Multiple Count columns"
            count_col = i
    assert len(list(filter(lambda x: x == -1, [mid_col, sub_col, \
        count_col]))) == 0, "missing column(s)"

def csv_aggregator():
    global matrix
    global line_count
    for row in csv_reader:
        if line_count > 0:
            aa1 = code[row[mid_col]]
            aa2 = code[row[sub_col]]
            count = int(row[count_col])
            if matrix[aa1][aa2] == None:
                matrix[aa1][aa2] = 1
            else:
                matrix[aa1][aa2] += count
            line_count += 1

# Writes the generated matrix to a CSV file
def write_to_csv(csv_file, bool):
    if not bool:
        return
    with open(csv_file, 'w') as f:
        top = ["", ""]
        writer = csv.writer(f)
        for i in range(20):
            top.append(list(code.keys())[i])
        writer.writerow(top)
        index = 0
        for line in matrix:
            writer.writerow(index, [list(code.keys())[index]] + [str(r) for r in line])
            index += 1
    print("CSV successfully written.")


csv_reader = csv.reader(open(sub_csv))
line_count = 0
mid_col = -1
sub_col = -1
count_col = -1
matrix = []

for i in range(20):
    lst = []
    for j in range(20):
        lst += [None]
    matrix.append(lst)

column_headers()
csv_aggregator()
write_to_csv(csv_file, generate_csv)

i = 0
for line in matrix:
    print(list(code.keys())[i] + str(line))
    i+=1