import itertools
import csv

'''making a master list of all possible combinations of 4 assays calls in ACM9 panel'''

#  use itertools to make every possible combination of 4 assays with homo/hemi/null placeholders
calls = ['Homo', 'Hemi', 'Null']
poss_combs = list(itertools.product(calls, repeat=4))

with open('tlm.csv', 'rb') as infile:
    reader = csv.reader(infile)

    with open('castle_master.csv', 'wb') as outfile:
        headers = ['DBSNP357223 Zygosity Call', 'DBSNP357253 Zygosity Call', 'DBSNP357255 Zygosity Call',
                   'DBSNP357256 Zygosity Call', 'ACM N9 Summary Confirm']
        writer = csv.DictWriter(outfile, fieldnames=headers)
        writer.writeheader()

        for line in reader:
            homo = 0
            hemi = 0
            null = 0
            no_call = 0
            no_data = 0

            for i in line:
                if i == 'Homo':
                    homo += 1
                if i == 'Hemi':
                    hemi += 1
                if i == 'Null':
                    null += 1
                if i == 'No Call':
                    no_call += 1
                if i == 'No Data':
                    no_data += 1

            if homo == 4:
                tlm_call = 'Trait'
            elif hemi == 4:
                tlm_call = 'Seg'
            elif null == 4:
                tlm_call = 'Wildtype'
            elif homo in range(1, 4) and hemi in range(1, 4) and null == 0:
                tlm_call = 'SegT'
            elif null in range(1, 4) and hemi in range(1, 4) and homo == 0:
                tlm_call = 'SegW'
            elif homo == 2 and null == 2:
                tlm_call = 'Undecided'
            elif homo == 3 and null == 1:
                tlm_call = 'BreakT'
            elif homo == 1 and null == 3:
                tlm_call = 'BreakW'
            else:
                tlm_call = 'Undecided'

            writer.writerow({'DBSNP357223 Zygosity Call': line[0], 'DBSNP357253 Zygosity Call': line[1],
                             'DBSNP357255 Zygosity Call': line[2], 'DBSNP357256 Zygosity Call': line[3],
                             'ACM N9 Summary Confirm': tlm_call})


