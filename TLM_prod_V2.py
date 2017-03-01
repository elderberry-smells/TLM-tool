import csv
import pandas as pd
from mas_class_tools import TraitN9
from mas_class_tools import TraitCastle
from mas_class_tools import TraitCR
import os

'''Version 2 of the TLM panel tool-- March 2017'''

kraken_file = 'SK_CrossingBlock PG 1419 RR B line and Crossing Block Set 1_no CR.csv'
find_the_dot = kraken_file.find('.')
final_name = kraken_file[:find_the_dot] + '_completed.xlsx'


# ---------------------------------------- TLM Panels ----------------------------------------------------------------

def determine_merge(list_panel):
    """determine the panel in the final merge sequence by assigning values of 1, 10, and 100 for N9, cas, and CR in
    the final_panel_list respectively.
    :param list_panel: is a list that has appended the TLM panels that exist in the file being analyzed
    Value determination is as follows:
    n9 = 1
    cas = 10
    cr = 100
    all three = 111
    n9 and cas = 11
    n9 and cr = 101
    cas and cr = 110
    """

    if 'n9' in list_panel:
        #  change the header to read Cyto Call
        list_panel = [1 if panel == 'n9' else panel for panel in list_panel]

    if 'castle' in list_panel:
        #  change the header to read GT200 Call
        list_panel = [10 if panel == 'castle' else panel for panel in list_panel]

    if 'cr' in list_panel:
        list_panel = [100 if panel == 'cr' else panel for panel in list_panel]

    panel_sum = 0
    for num in list_panel:
        panel_sum += num

    return panel_sum


#  Determine what assays were run in the Kraken study, and if anything needs to be analyzed of converted in the file

with open(kraken_file, 'rb') as kraken:
    kreader = csv.DictReader(kraken)
    kheaders = kreader.fieldnames

kdf = pd.read_csv(kraken_file)

pandas_list = []  # an empty list to append each panels dataframe name for a final merge sequence

# ----------------------------------------- N9 Panel ------------------------------------------------------------------

"""determine if N9 panel is in the kraken study, and if so makes a df with a N9 call"""

n9_panel = ['DBSNP357223 Zygosity Call', 'DBSNP357253 Zygosity Call', 'DBSNP357255 Zygosity Call',
            'DBSNP357256 Zygosity Call']

#  index the headers to see if the N9 panel is in the kraken file, and if so, which column #'s they are in
n9_index = []
for i in n9_panel:
    try:
        n9_index.append(kheaders.index(i))
    except:
        continue
n9_len = int(len(n9_index))

# if N9 panel is in the kraken file, make a pandas dataframe, and write that section to a new csv to run analysis
if n9_len > 0:
    n9_headers = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']
    for i in n9_index:
        n9_headers.append(kheaders[i])
    n9df = kdf[n9_headers]
    n9df.to_csv('n9_temp.csv', index=False)

    with open('n9_temp.csv', 'rb') as n9input:
        n9reader = csv.DictReader(n9input)
        n9headers = n9reader.fieldnames
        n9headers.append('ACM N9 Summary Confirm')

        with open('n9_call.csv', 'wb') as n9output:
            n9writer = csv.DictWriter(n9output, fieldnames=n9headers)
            n9writer.writeheader()

            for line in n9reader:
                n9id = TraitN9(line['DBSNP357223 Zygosity Call'], line['DBSNP357253 Zygosity Call'],
                               line['DBSNP357255 Zygosity Call'], line['DBSNP357256 Zygosity Call'])
                n9_call = n9id.n9_call()

                n9writer.writerow({'Box': line['Box'], 'Well': line['Well'], 'Project': line['Project'],
                                   'Pedigree': line['Pedigree'], 'Source ID': line['Source ID'],
                                   'Geno_Id': line['Geno_Id'], 'RowId': line['RowId'], 'Loc Seq#': line['Loc Seq#'],
                                   'DBSNP357223 Zygosity Call': line['DBSNP357223 Zygosity Call'],
                                   'DBSNP357253 Zygosity Call': line['DBSNP357253 Zygosity Call'],
                                   'DBSNP357255 Zygosity Call': line['DBSNP357255 Zygosity Call'],
                                   'DBSNP357256 Zygosity Call': line['DBSNP357256 Zygosity Call'],
                                   'ACM N9 Summary Confirm': n9_call})

            pandas_list.append('n9')

    os.remove('n9_temp.csv')

# --------------------------------------- Castle Panel ----------------------------------------------------------------

"""determine if Castle panel is in the kraken study, and if so makes a df with a Castle call"""

castle_panel = ['298512 Zygosity Call', '298516 Zygosity Call', '298518 Zygosity Call', '298535 Zygosity Call']

#  index the headers to see if the castle panel is in the kraken file, and if so, which column #'s they are in
castle_index = []
for i in castle_panel:
    try:
        castle_index.append(kheaders.index(i))
    except:
        continue
cas_len = int(len(castle_index))

# if that index is not empty, make a pandas dataframe, and write that section to a new csv to run analysis
if cas_len > 0:
    castle_headers = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']
    for i in castle_index:
        castle_headers.append(kheaders[i])
    casdf = kdf[castle_headers]
    casdf.to_csv('cas_temp.csv', index=False)

    with open('cas_temp.csv', 'rb') as casinput:
        casreader = csv.DictReader(casinput)
        casheaders = casreader.fieldnames
        casheaders.append('Castle Summary')

        with open('cas_call.csv', 'wb') as casoutput:
            caswriter = csv.DictWriter(casoutput, fieldnames=casheaders)
            caswriter.writeheader()

            for line in casreader:
                casid = TraitCastle(line['298512 Zygosity Call'], line['298516 Zygosity Call'],
                                    line['298518 Zygosity Call'], line['298535 Zygosity Call'])
                cas_call = casid.cas_call()

                caswriter.writerow({'Box': line['Box'], 'Well': line['Well'], 'Project': line['Project'],
                                    'Pedigree': line['Pedigree'], 'Source ID': line['Source ID'],
                                    'Geno_Id': line['Geno_Id'], 'RowId': line['RowId'], 'Loc Seq#': line['Loc Seq#'],
                                    '298512 Zygosity Call': line['298512 Zygosity Call'],
                                    '298516 Zygosity Call': line['298516 Zygosity Call'],
                                    '298518 Zygosity Call': line['298518 Zygosity Call'],
                                    '298535 Zygosity Call': line['298535 Zygosity Call'],
                                    'Castle Summary': cas_call})

            pandas_list.append('castle')

    os.remove('cas_temp.csv')

# ----------------------------------------- CR Panel ------------------------------------------------------------------
"""determine if Castle panel is in the kraken study, and if so makes a df with a Castle call"""

cr_panel = ['CRM1 Zygosity Call', 'CRM2 Zygosity Call']

#  index the headers to see if the CR panel is in the kraken file, and if so, which column #'s they are in
cr_index = []
for i in cr_panel:
    try:
        cr_index.append(kheaders.index(i))
    except:
        continue
cr_len = int(len(cr_index))

# if that index is not empty, make a pandas dataframe, and write that section to a new csv to run analysis
if cr_len > 0:
    cr_headers = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']
    for i in cr_index:
        cr_headers.append(kheaders[i])
    crdf = kdf[cr_headers]
    crdf.to_csv('cr_temp.csv', index=False)

    with open('cr_temp.csv', 'rb') as crinput:
        crreader = csv.DictReader(crinput)
        crheaders = crreader.fieldnames
        crheaders.append('CR Mendel Summary Confirm')

        with open('cr_call.csv', 'wb') as croutput:
            crwriter = csv.DictWriter(croutput, fieldnames=crheaders)
            crwriter.writeheader()

            for line in crreader:
                crid = TraitCR(line['CRM1 Zygosity Call'], line['CRM2 Zygosity Call'])
                cr_call = crid.cr_call()

                crwriter.writerow({'Box': line['Box'], 'Well': line['Well'], 'Project': line['Project'],
                                   'Pedigree': line['Pedigree'], 'Source ID': line['Source ID'],
                                   'Geno_Id': line['Geno_Id'], 'RowId': line['RowId'], 'Loc Seq#': line['Loc Seq#'],
                                   'CRM1 Zygosity Call': line['CRM1 Zygosity Call'],
                                   'CRM2 Zygosity Call': line['CRM2 Zygosity Call'],
                                   'CR Mendel Summary Confirm': cr_call})

            pandas_list.append('cr')

    os.remove('cr_temp.csv')

# ----------------------------------------- Merge with Pandas --------------------------------------------------------

# read each panels analyzed data as a dataframe, and merge them in sequence.

krak_info = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']

df_final = pd.read_csv(kraken_file)  # Kraken dataframe
df_initial = df_final[krak_info]  # sample and well info (no molecular data)

merge_number = determine_merge(pandas_list)

#  begin merging in sequence.  start with initial df, merging the available df's, and end with merge on the final df

if merge_number == 1:  # if the panel consists of only N9
    n9_df = pd.read_csv('n9_call.csv')

    merge_init = pd.merge(left=df_initial, right=n9_df, how='left')
    merge_final = pd.merge(left=merge_init, right=df_final, how='left')
    merge_final.to_excel(final_name, index=False)
    os.remove('n9_call.csv')

elif merge_number == 10:  # if the panel consists of only Castle
    cas_df = pd.read_csv('cas_call.csv')

    merge_init = pd.merge(left=df_initial, right=cas_df, how='left')
    merge_final = pd.merge(left=merge_init, right=df_final, how='left')
    merge_final.to_excel(final_name, index=False)
    os.remove('cas_call.csv')


elif merge_number == 11:  # if the panel consists of N9 and Castle
    n9_df = pd.read_csv('n9_call.csv')
    cas_df = pd.read_csv('cas_call.csv')

    merge_init = pd.merge(left=df_initial, right=n9_df, how='left')
    merge_two = pd.merge(left=merge_init, right=cas_df, how='left')
    merge_final = pd.merge(left=merge_two, right=df_final, how='left')
    merge_final.to_excel(final_name, index=False)
    os.remove('n9_call.csv')
    os.remove('cas_call.csv')

elif merge_number == 100:  # if the panel consists of only CR
    cr_df = pd.read_csv('cr_call.csv')

    merge_init = pd.merge(left=df_initial, right=cr_df, how='left')
    merge_final = pd.merge(left=merge_init, right=df_final, how='left')
    merge_final.to_excel(final_name, index=False)
    os.remove('cr_call.csv')

elif merge_number == 101:  # if the panel consists of N9 and CR
    n9_df = pd.read_csv('n9_call.csv')
    cr_df = pd.read_csv('cr_call.csv')

    merge_init = pd.merge(left=df_initial, right=n9_df, how='left')
    merge_two = pd.merge(left=merge_init, right=cr_df, how='left')
    merge_final = pd.merge(left=merge_two, right=df_final, how='left')
    merge_final.to_excel(final_name, index=False)
    os.remove('n9_call.csv')
    os.remove('cr_call.csv')

elif merge_number == 110:  # if the panel consists of Castle and CR
    cas_df = pd.read_csv('cas_call.csv')
    cr_df = pd.read_csv('cr_call.csv')

    merge_init = pd.merge(left=df_initial, right=cas_df, how='left')
    merge_two = pd.merge(left=merge_init, right=cr_df, how='left')
    merge_final = pd.merge(left=merge_two, right=df_final, how='left')
    merge_final.to_excel(final_name, index=False)
    os.remove('cas_call.csv')
    os.remove('cr_call.csv')

elif merge_number == 111:  # if the panel consists of N9, Castle, and CR
    n9_df = pd.read_csv('n9_call.csv')
    cas_df = pd.read_csv('cas_call.csv')
    cr_df = pd.read_csv('cr_call.csv')

    merge_init = pd.merge(left=df_initial, right=n9_df, how='left')
    merge_two = pd.merge(left=merge_init, right=cas_df, how='left')
    merge_three = pd.merge(left=merge_two, right=cr_df, how='left')
    merge_final = pd.merge(left=merge_three, right=df_final, how='left')
    merge_final.to_excel(final_name, index=False)
    os.remove('n9_call.csv')
    os.remove('cas_call.csv')
    os.remove('cr_call.csv')

else:
    merge_final = pd.merge(left=df_initial, right=df_final, how='left')
    merge_final.to_excel(final_name, index=False)
