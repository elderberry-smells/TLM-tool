import itertools
import csv

possible = ['Homo', 'Hemi', 'Null', 'No Call', 'No Data']

other_business = list(itertools.product(possible, repeat = 2))

with open('cr_master.csv', 'wb') as outfile:
    headers = ['CRM1 Zygosity Call','CRM2 Zygosity Call', 'CR Mendel Summary Confirm']
    writer = csv.DictWriter(outfile, fieldnames=headers)
    writer.writeheader()

    for line in other_business:
        writer.writerow({'CRM1 Zygosity Call':line[0], 'CRM2 Zygosity Call': line[1],
                         'CR Mendel Summary Confirm': 'confirm'})



