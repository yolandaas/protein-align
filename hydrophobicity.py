from math import *
import csv
from aminoacids import *

###############
# DESCRIPTION #
###############
"""
This file serves to generate 3 matrices based on substitution 
frequencies, given hydrophobic/philic flanks. The script requires 
import from a substitution matrix csv file generated by matrixgen.py

Author Yolanda Shen

IN PROGRESS
"""

# Variables to be modified
sub_csv = 'C:/Users/yolan/SABmark/SABmark/sup/sub_matrix.csv'
hydrophob_csv_file = 'C:/Users/yolan/SABmark/SABmark/sup/hydrophob_matrix.csv'
hydrophil_csv_file = 'C:/Users/yolan/SABmark/SABmark/sup/hydrophil_matrix.csv'
combo_csv_file = 'C:/Users/yolan/SABmark/SABmark/sup/combo_hydro_matrix.csv'
generate_hydrophob_csv = False
generate_hydrophil_csv = False
generate_combo_csv = False

#############################################
# Do not modify anything below this comment #
#############################################

def csv_aggregator():
    global hydrophob_matrix
    csv_reader.__next__()
    for row in csv_reader:
        lflank = row[0]
        mid = row[1]
        rflank = row[2]
        sub = row[3]
        count = int(row[4])
        key = mid + sub
        if lflank in hydrophob and rflank in hydrophob:
            if key not in hydrophob_matrix.keys() :
                hydrophob_matrix[key] = count
            else:
                hydrophob_matrix[key] += count
        elif lflank not in hydrophob and rflank not in hydrophob:
            if key not in hydrophil_matrix.keys() :
                hydrophil_matrix[key] = count
            else:
                hydrophil_matrix[key] += count
        else:
            if key not in combo_matrix.keys() :
                combo_matrix[key] = count
            else:
                combo_matrix[key] += count

# Writes the generated matrix to a CSV file
def write_to_csv(matrix, csv_file, bool):
    if not bool:
        return
    with open(csv_file, 'w') as f:
        for key in matrix.keys():
            if key == 'type':
                f.write("Parent,Sub,Count")
            else:
                f.write("%s,%s,%s,\n" % (key[0], key[1], matrix[key]))
    print("CSV successfully written.")


csv_reader = csv.reader(open(sub_csv))
hydrophob_matrix = {'type':'hydrophob'}
hydrophil_matrix = {'type':'hydrophil'}
combo_matrix = {'type':'combo'}

csv_aggregator()
write_to_csv(hydrophob_matrix, hydrophob_csv_file, generate_hydrophob_csv)
write_to_csv(hydrophil_matrix, hydrophil_csv_file, generate_hydrophil_csv)
write_to_csv(combo_matrix, combo_csv_file, generate_combo_csv)

print(hydrophob_matrix)

