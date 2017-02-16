import csv
import pandas as pd
import os

'''Version 1 of the TLM panel tool-- February 2017'''

kraken_file = 'SK_CrossingBlock PG 1419 RR B line and Crossing Block Set 1.csv'
find_the_dot = kraken_file.find('.')
final_name = kraken_file[:find_the_dot] + '_completed.xlsx'

#  get a dictionary of each master file for the TLM panels (N9, Castle, CR)

with open('n9_master.csv') as master:
    n9reader = csv.reader(master)
    n9id = dict((rows[4], rows[0:4]) for rows in n9reader)

with open('castle_master.csv') as master:
    casreader = csv.reader(master)
    casid = dict((rows[4], rows[0:4]) for rows in casreader)

with open('cr_master.csv') as master:
    crreader = csv.reader(master)
    crid = dict((rows[2], rows[0:2]) for rows in crreader)

#  identify the panel of assays used for each TLM, and their associated headers from a kraken file

n9_panel = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#',
            'DBSNP357223 Zygosity Call', 'DBSNP357253 Zygosity Call', 'DBSNP357255 Zygosity Call',
            'DBSNP357256 Zygosity Call']

castle_panel = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#',
                '298512 Zygosity Call', '298516 Zygosity Call', '298518 Zygosity Call', '298535 Zygosity Call']

cr_panel = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#',
            'CRM1 Zygosity Call', 'CRM2 Zygosity Call']

#  make a dataframe for each panel and save it as a temp file.  Run the profile script on that temp file.

df = pd.read_csv(kraken_file)  # Kraken dataframe
n9df = df[n9_panel]  # N9 dataframe
casdf = df[castle_panel]  # castle dataframe
crdf = df[cr_panel]  # clubroot dataframe

n9df.to_csv('n9_temp.csv', index=False)
casdf.to_csv('cas_temp.csv', index=False)
crdf.to_csv('cr_temp.csv', index=False)

# open the N9 temporary file and add the call referencing the master id, save as a new "call file", delete the temp file
with open('n9_temp.csv', 'r') as infile:
    dreader = csv.DictReader(infile)
    headers = dreader.fieldnames
    headers.append('ACM N9 Summary Confirm')

    with open('n9_call.csv', 'wb') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=headers)
        writer.writeheader()

        for line in dreader:
            n9_assid = [line['DBSNP357223 Zygosity Call'], line['DBSNP357253 Zygosity Call'],
                        line['DBSNP357255 Zygosity Call'], line['DBSNP357256 Zygosity Call']]

            try:
                keys = n9id.keys()[n9id.values().index(n9_assid)]
                end_id = keys.find("(")
                profile_id = keys[:end_id]
            except:
                profile_id = 'No Call'

            writer.writerow({
                'Box': line['Box'], 'Well': line['Well'], 'Project': line['Project'],
                'Pedigree': line['Pedigree'], 'Source ID': line['Source ID'], 'Geno_Id': line['Geno_Id'],
                'RowId': line['RowId'], 'Loc Seq#': line['Loc Seq#'],
                'DBSNP357223 Zygosity Call': line['DBSNP357223 Zygosity Call'],
                'DBSNP357253 Zygosity Call': line['DBSNP357253 Zygosity Call'],
                'DBSNP357255 Zygosity Call': line['DBSNP357255 Zygosity Call'],
                'DBSNP357256 Zygosity Call': line['DBSNP357256 Zygosity Call'],
                "ACM N9 Summary Confirm": profile_id})
os.remove('n9_temp.csv')

# open the Castle temporary file and add the call referencing the master id, save as a new "call file",
# delete the temp file
with open('cas_temp.csv', 'r') as infile:
    casreader = csv.DictReader(infile)
    casheaders = casreader.fieldnames
    casheaders.append('Castle Summary')

    with open('cas_call.csv', 'wb') as outfile:
        caswriter = csv.DictWriter(outfile, fieldnames=casheaders)
        caswriter.writeheader()

        for line in casreader:
            cas_assid = [line['298512 Zygosity Call'], line['298516 Zygosity Call'],
                         line['298518 Zygosity Call'], line['298535 Zygosity Call']]

            try:
                keys = casid.keys()[casid.values().index(cas_assid)]
                end_id = keys.find("(")
                cas_id = keys[:end_id]
            except:
                cas_id = 'No Call'

            caswriter.writerow({
                'Box': line['Box'], 'Well': line['Well'], 'Project': line['Project'],
                'Pedigree': line['Pedigree'], 'Source ID': line['Source ID'], 'Geno_Id': line['Geno_Id'],
                'RowId': line['RowId'], 'Loc Seq#': line['Loc Seq#'],
                '298512 Zygosity Call': line['298512 Zygosity Call'],
                '298516 Zygosity Call': line['298516 Zygosity Call'],
                '298518 Zygosity Call': line['298518 Zygosity Call'],
                '298535 Zygosity Call': line['298535 Zygosity Call'],
                'Castle Summary': cas_id})
os.remove('cas_temp.csv')

# open the CR temporary file and add the call referencing the master id, save as a new "call file", delete the temp file
with open('cr_temp.csv', 'r') as crfile:
    crreader = csv.DictReader(crfile)
    crheaders = crreader.fieldnames
    crheaders.append('CR Mendel Summary Confirm')

    with open('cr_call.csv', 'wb') as croutfile:
        crwriter = csv.DictWriter(croutfile, fieldnames=crheaders)
        crwriter.writeheader()

        for line in crreader:
            cr_assid = [line['CRM1 Zygosity Call'], line['CRM2 Zygosity Call']]

            try:
                keys = crid.keys()[crid.values().index(cr_assid)]
                end_id = keys.find("(")
                cr_id = keys[:end_id]
            except:
                cr_id = 'No Call'

            crwriter.writerow({
                'Box': line['Box'], 'Well': line['Well'], 'Project': line['Project'],
                'Pedigree': line['Pedigree'], 'Source ID': line['Source ID'], 'Geno_Id': line['Geno_Id'],
                'RowId': line['RowId'], 'Loc Seq#': line['Loc Seq#'],
                'CRM1 Zygosity Call': line['CRM1 Zygosity Call'],
                'CRM2 Zygosity Call': line['CRM2 Zygosity Call'],
                'CR Mendel Summary Confirm': cr_id})

os.remove('cr_temp.csv')

# Read all the call files, and merge them into each other, also merge the original Kraken file
df_n9 = pd.read_csv('n9_call.csv')  # N9 dataframe with call
df_cas = pd.read_csv('cas_call.csv')  # Castle dataframe with call
df_cr = pd.read_csv('cr_call.csv')  # Clubroot dataframe with call

# merge N9 and Castle
merge_one = pd.merge(left=df_n9, right=df_cas, how='left')

# merge merge_one and CR
merge_two = pd.merge(left=merge_one, right=df_cr, how='left')

# merge merge_two and original dataframe (kraken file) and write it as final csv
merge_final = pd.merge(left=merge_two, right=df, how='left')

# remove temp files
os.remove('n9_call.csv')
os.remove('cas_call.csv')
os.remove('cr_call.csv')

# Write the final merged file to an excel, naming it as the original + "_completed"
merge_final.to_excel(final_name, index=False)
