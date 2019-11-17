from pathlib import Path
from access import *

# Main area to run tests
path_to_sup_test = 'C:/Users/yolan/SABmark/SABmark/sup/'
iter_test = False
nseq_test = True
run_in_test = False
run_del_test = False
run_sub_test = False

# file_iterator test
if iter_test:
    group = 1
    test = glob.iglob(path_to_sup_test + 'group' + str(group) + '/reference/*.fasta')
    x = file_iterator("sub", test, group)
    print(matrix)
    print(x)
    test1 = glob.iglob(path_to_sup_test + 'group1/reference/*.fasta')
    file_iterator("del", test1)


# nseq test 
nseq = gen_nseq()
for g in nseq.keys():
    if not nseq_test:
        break
    print(str(g) + " " + str(nseq[g]))


# neighbor test
x = "VLS----------------------------EGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKAS-----------EDLKKHGVTVLTALGAILKKKGH----HEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGY-"
y = "---SIVTKSIVNADAEARYLSPGELDRIKSFVTSGERRVRIAETMTGARERIIKQAGDQLFGKRPDVVS-----------------PGGNAYGADMTATCLRDLDYYLRLITYGIVAGD-VTPIEEIGVVGVREMYKSLGTP-IEAIAEGVRAMKSVATSL----LSGADAAEAGSYFDYLIGAMS---------"
assert(neighbor(x, 0) == ['$', 'V', 'L', 1])
assert(neighbor(x, 2) == ['L', 'S', 'E', 29])
assert(neighbor(y, 2) == ['$', '-', 'S', 3])
assert(neighbor(x, 3) == ['S', '-', 'E', 28,])
assert(neighbor(y, 3) == ['$', 'S', 'I', 4])

# in_matrix test
if run_in_test:
    print(in_matrix(['B---Z', 'AHDJA']))
    print(in_matrix(['----Z', 'AHDJA']))
    print(in_matrix(['B----', 'AHDJA']))
    print(in_matrix([x, y]))

# del_matrix test
if run_del_test:
    print(del_matrix(['AHDJA', 'B---Z']))
    print(del_matrix(['AHDJA', '----Z']))
    print(del_matrix(['AHDJA', 'B----']))
    print(del_matrix([x, y]))

# sub_matrix test
if run_sub_test:
    print(sub_matrix(['AAAAAA', 'DADBAC'])) 
        # should be {'type': 'sub', '$AAD': 1, 'AAAD': 1, 'AAAB': 1, 'AA$C': 1}, without the 'type'

    test_sub = glob.iglob(path_to_sup_test + '/group1/reference/*.fasta').__iter__()
    print(sub_matrix(get_seqs(test_sub.__next__())))