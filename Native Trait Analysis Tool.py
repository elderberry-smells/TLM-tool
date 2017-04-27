import csv
import pandas as pd
from pandas import ExcelWriter
from TLM import TraitLinkedMarker as Trait
from summary_class import *
import glob
import os
import time

'''Version 3 of the TLM panel tool-- Validated April 2017 -- [Author: Brian James]'''

#  get a list of all the csv files to run the script on
path = r'C:\Users\U590135\Desktop\Python\tlm_tool\TLM\Version 3'
extension = 'csv'
os.chdir(path)
result = [i for i in glob.glob('*.{}'.format(extension))]

completed_path = r'C:\Users\U590135\Desktop\Python\tlm_tool\TLM\Version 3\completed'

results_files = []
for i in result:
    if 'temp_' in i:
        continue
    if '_call' in i:
        continue
    else:
        results_files.append(i)

#  ----------------------------------  read the folder and find files for analysis  -----------------------------------

for fname in results_files:
    print('Analyzing ' + fname + '...')

    #  make the naming conventions for the files being written in the end
    encodename = 'encoded_' + fname
    find_the_dot = fname.find('.')
    final_name = fname[:find_the_dot] + '_Report.xlsx'
    results_writer = ExcelWriter(final_name, options={'encoding': 'utf-8'})

    #  ---------------------------------------  Encoding the file to UTF-8  -------------------------------------------
    with open(fname, 'r+') as encode_kraken:
        encode_reader = csv.DictReader(encode_kraken)
        encode_headers = encode_reader.fieldnames

    with open(encodename, 'wb') as encode_file:  # write a temporary UTF-8 encoded version of the kraken file
        writer = csv.DictWriter(encode_file, fieldnames=encode_headers)

        r = open(fname, 'rb')
        remove_controls = ['H7', 'H07', 'H8', 'H08', 'H9', 'H09', 'H10', 'H11', 'H12']

        for line in r:
            bits = line.split(',')
            if bits[1] in remove_controls:
                continue
            else:
                bits[3] = bits[3].decode('latin-1').encode('utf-8')  # encode the pedigree col as UTF-8, not Latin-1
                encode_file.write(','.join(bits))

        r.close()

    #  ----------------------------------------  Summary Table  --------------------------------------------------------

    with open(encodename, 'rb') as summary_kraken:
        sum_reader = csv.DictReader(summary_kraken)
        sum_headers = sum_reader.fieldnames
        sum_assays = [i for i in sum_headers if 'Zygosity' in i]  # only assays in the file for the summary

        summary_init = SummaryTable(sum_reader, sum_headers)  # initiate the class method
        summary_dict = summary_init.get_summary()  # get the summary dictionary from the SummaryTable class
        summary_df1 = pd.DataFrame.from_dict(summary_dict)  # convert the dictionary to a dataframe in pandas
        summary_df2 = summary_df1[sum_assays]  # make sure the columns match the original kraken file
        summary_df = summary_df2.reindex(["Trait", "Seg", "Wildtype",
                                          "No Call", "Fail", "% Data Return"])  # order the index

    # ----------------------------------- Analysis of the Kraken File -------------------------------------------------

    #  Determine the assays run in the Kraken study, and if anything needs to be analyzed or converted in the file
    with open(encodename, 'rb') as kraken:
        kreader = csv.DictReader(kraken)
        kheaders = kreader.fieldnames
        creader = csv.reader(kraken)

        #  ---------------------------- Cyto, GT200, BGSP-02 and LepR3 Conversion ------------------------------------

        convert_df = pd.read_csv(encodename)

        if 'Cyto Zygosity Call' in convert_df.columns:
            convert_df = convert_df.replace(to_replace={'Cyto Zygosity Call': {'Trait': 'A-Cyto'}})
            convert_df = convert_df.replace(to_replace={'Cyto Zygosity Call': {'Wildtype': 'B-Cyto'}})
            convert_df.rename(columns={'Cyto Zygosity Call': 'Cyto Call'}, inplace=True)

        if 'GT 200 Zygosity Call' in convert_df.columns:
            convert_df = convert_df.replace(to_replace={'GT 200 Zygosity Call': {'Homo': 'Plus'}})
            convert_df = convert_df.replace(to_replace={'GT 200 Zygosity Call': {'Null': 'Minus'}})
            convert_df.rename(columns={'GT 200 Zygosity Call': 'GT200 Call'}, inplace=True)

        if 'GT200 Zygosity Call' in convert_df.columns:
            convert_df = convert_df.replace(to_replace={'GT200 Zygosity Call': {'Homo': 'Plus'}})
            convert_df = convert_df.replace(to_replace={'GT200 Zygosity Call': {'Null': 'Minus'}})
            convert_df.rename(columns={'GT200 Zygosity Call': 'GT200 Call'}, inplace=True)

        if 'LepR3-R3S3A Zygosity Call' in convert_df.columns:
            convert_df = convert_df.replace(to_replace={'LepR3-R3S3A Zygosity Call': {'Homo': 'Trait'}})
            convert_df = convert_df.replace(to_replace={'LepR3-R3S3A Zygosity Call': {'Hemi': 'Seg'}})
            convert_df = convert_df.replace(to_replace={'LepR3-R3S3A Zygosity Call': {'Null': 'Wildtype'}})

        if 'BGSP-02 Zygosity Call' in convert_df.columns:
            convert_df = convert_df.replace(to_replace={'BGSP-02 Zygosity Call': {'Trait': 'Plus'}})
            convert_df = convert_df.replace(to_replace={'BGSP-02 Zygosity Call': {'Wildtype': 'Minus'}})
            convert_df.rename(columns={'BGSP-02 Zygosity Call': 'BGSP-02'}, inplace=True)

        kdf = convert_df
        df_final = convert_df

    pandas_list = []  # an empty list to append each panels dataframe name (if used) for a final merge sequence

    # ----------------------------------------- N9 Panel --------------------------------------------------------------

    """determine if N9 panel is in the kraken study, and if so make a dataframe with a N9 call"""

    n9_panel = ['DBSNP357223 Zygosity Call', 'DBSNP357226 Zygosity Call', 'DBSNP357233 Zygosity Call',
                'DBSNP357252 Zygosity Call', 'DBSNP357253 Zygosity Call', 'DBSNP357255 Zygosity Call',
                'DBSNP357256 Zygosity Call', 'DBSNP357257 Zygosity Call']

    #  index the headers to see if the N9 panel is in the kraken file, and if so, which column #'s they are in
    n9_index = []
    for i in n9_panel:
        # noinspection PyBroadException
        try:
            n9_index.append(kheaders.index(i))
        except:
            continue
    n9_len = int(len(n9_index))

    # if N9 panel is in the kraken file, make a dataframe, and write that section to a new csv to run analysis
    if n9_len >= 2:
        #  keep the unique identifiers from kraken study in new temp dataframe
        n9_headers = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']

        for i in n9_index:  # add the N9 assays in the kraken file to the header list
            n9_headers.append(kheaders[i])

        n9df = kdf[n9_headers]  # make a dataframe of the N9 panel as a temp file
        n9df.to_csv('n9_temp.csv', index=False)

        with open('n9_temp.csv', 'rb') as n9input:  # open newly created dataframe of N9 panel
            n9reader = csv.DictReader(n9input)
            n9headers = n9reader.fieldnames

            with open('n9_call.csv', 'wb') as n9output:  # write a new call file with the N9 call column added
                n9headers.append('ACM N9 Summary')
                n9writer = csv.DictWriter(n9output, fieldnames=n9headers)
                n9writer.writeheader()

                for line in n9reader:
                    dict_n9 = {}  # make a dictionary for each line in temp file, to pass into the TLM class
                    if len(n9_index) > 0:
                        for i in n9_index:
                            dict_n9[kheaders[i]] = line[kheaders[i]]

                    n9 = Trait(dict_n9)  # initiate the TLM class
                    n9_call = n9.get_tlm_call()  # get the trait/seg/wildtype call for this line in the file
                    line['ACM N9 Summary'] = n9_call  # add call to the DictReader line

                    n9writer.writerow(line)  # write the line (with the call now in the line) to the temp call file

                pandas_list.append('n9')  # add n9 (string) to pandas list for prompts to merge later

        os.remove('n9_temp.csv')

    # --------------------------------------- Castle Panel ------------------------------------------------------------

    """determine if Castle panel is in the kraken study, and if so makes a df with a Castle call"""

    castle_panel = ['298512 Zygosity Call', '298516 Zygosity Call', '298518 Zygosity Call', '298535 Zygosity Call']

    #  index the headers to see if the castle panel is in the kraken file, and if so, which column #'s they are in
    castle_index = []
    for i in castle_panel:
        # noinspection PyBroadException
        try:
            castle_index.append(kheaders.index(i))
        except:
            continue
    cas_len = int(len(castle_index))

    # if that index is not empty, make a pandas dataframe, and write that section to a new csv to run analysis
    if cas_len >= 2:
        #  keep the unique identifiers from kraken study in new temp dataframe
        castle_headers = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']

        for i in castle_index:  # add the Castle assays in the kraken file to the header list
            castle_headers.append(kheaders[i])

        casdf = kdf[castle_headers]
        casdf.to_csv('cas_temp.csv', index=False)  # make a temp dataframe for just castle data

        with open('cas_temp.csv', 'rb') as casinput:  # read the newly created dataframe
            csreader = csv.DictReader(casinput)
            csheaders = csreader.fieldnames

            with open('cas_call.csv', 'wb') as casoutput:  # create another temporary file for the TLM call
                csheaders.append('Castle Summary')  # add a column for the castle call
                cswriter = csv.DictWriter(casoutput, fieldnames=csheaders)
                cswriter.writeheader()

                for line in csreader:
                    dict_cas = {}  # make a dictionary for each line in temp file, to pass into the TLM class
                    if len(castle_index) > 0:
                        for i in castle_index:
                            dict_cas[kheaders[i]] = line[kheaders[i]]

                    castle = Trait(dict_cas)  # initiate the TLM class with the dictionary
                    castle_call = castle.get_tlm_call()  # get the call for each line
                    line['Castle Summary'] = castle_call  # add call to the DictReader line

                    cswriter.writerow(line)  # write the line (with the call now in the line) to the temp call file

                pandas_list.append('castle')  # add castle (string) to pandas list for prompts to merge later

        os.remove('cas_temp.csv')

    # ----------------------------------------- CR Panel --------------------------------------------------------------
    """determine if CR panel is in the kraken study, and if so makes a df with a CR call"""

    cr_panel = ['CRM1 Zygosity Call', 'CRM2 Zygosity Call']

    #  index the headers to see if the CR panel is in the kraken file, and if so, which column #'s they are in
    cr_index = []
    for i in cr_panel:
        # noinspection PyBroadException
        try:
            cr_index.append(kheaders.index(i))
        except:
            continue
    cr_len = int(len(cr_index))

    # if that index is not empty, make a pandas dataframe, and write that section to a new csv to run analysis
    if cr_len >= 2:
        #  keep the unique identifiers from kraken study in new temp dataframe
        cr_headers = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']

        for i in cr_index:  # add the CR assays in the kraken file to the header list
            cr_headers.append(kheaders[i])

        crdf = kdf[cr_headers]  # make a dataframe of just CR data
        crdf.to_csv('cr_temp.csv', index=False)

        with open('cr_temp.csv', 'rb') as crinput:  # open that newly created dataframe
            crreader = csv.DictReader(crinput)
            crheaders = crreader.fieldnames

            with open('cr_call.csv', 'wb') as croutput:  # make a new file with the CR call column added
                crheaders.append('CR Mendel Summary')
                crwriter = csv.DictWriter(croutput, fieldnames=crheaders)
                crwriter.writeheader()

                for line in crreader:

                    dict_cr = {}  # make a dictionary for each line in temp file, to pass into the TLM class
                    if len(cr_index) > 0:
                        for i in cr_index:
                            dict_cr[kheaders[i]] = line[kheaders[i]]

                    cr = Trait(dict_cr)  # Initiate the TLM class with dictionary
                    cr_call = cr.get_tlm_call()  # get the trait/seg/wildtype call for this line in the file
                    line['CR Mendel Summary'] = cr_call  # add call to the DictReader line

                    crwriter.writerow(line)  # write the line (with the call now in the line) to the temp call file

                pandas_list.append('cr')  # add cr (string) to pandas list for prompts to merge later

        os.remove('cr_temp.csv')

    # -------------------------------------------  RLM7 Spring Analysis -----------------------------------------------
    """determine if RLM7 Spring panel is in the kraken study, and if so makes a df with a RLM7 spring call"""

    rlm7s_panel = ['295698 Zygosity Call', '295699 Zygosity Call']

    #  index the headers to see if the CR panel is in the kraken file, and if so, which column #'s they are in
    rlm7s_index = []
    for rlm7si in rlm7s_panel:
        # noinspection PyBroadException
        try:
            rlm7s_index.append(kheaders.index(rlm7si))
        except:
            continue
    rlm7s_len = int(len(rlm7s_index))

    # if that index is not empty, make a pandas dataframe, and write that section to a new csv to run analysis
    if rlm7s_len >= 2:
        #  keep the unique identifiers from kraken study in new temp dataframe
        rlm7s_headers = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']

        for i in rlm7s_index:  # add the RLM7 spring assays in the kraken file to the header list
            rlm7s_headers.append(kheaders[i])

        rlm7sdf = kdf[rlm7s_headers]  # make a dataframe of just RLM7 spring data
        rlm7sdf.to_csv('rlm7s_temp.csv', index=False)

        with open('rlm7s_temp.csv', 'rb') as rlm7sinput:  # open that newly created dataframe
            rlm7sreader = csv.DictReader(rlm7sinput)
            rlm7sheaders = rlm7sreader.fieldnames

            with open('rlm7s_call.csv', 'wb') as rlm7soutput:    # make a new file with the rlm7 call column added
                rlm7sheaders.append('RLM 7 Spring Summary')
                rlm7swriter = csv.DictWriter(rlm7soutput, fieldnames=rlm7sheaders)
                rlm7swriter.writeheader()

                for line in rlm7sreader:

                    dict_rlm7s = {}  # make a dictionary for each line in temp file, to pass into the TLM class
                    if len(rlm7s_index) > 0:
                        for i in rlm7s_index:
                            dict_rlm7s[kheaders[i]] = line[kheaders[i]]

                    rlm7s = Trait(dict_rlm7s)  # initiate the TLM class with dictionary
                    rlm7s_call = rlm7s.get_tlm_call()  # get the trait/seg/wildtype call for this line in the file
                    line['RLM 7 Spring Summary'] = rlm7s_call  # add call to the DictReader line

                    rlm7swriter.writerow(line)  # write the line (with the call now in the line) to the temp call file

                pandas_list.append('rlm7s')  # add rlm7s (string) to pandas list for prompts to merge later

        os.remove('rlm7s_temp.csv')

    # -------------------------------------------  RLM7 Winter Analysis -----------------------------------------------
    """determine if RLM7 Winter panel is in the kraken study, and if so makes a df with a RLM7 winter call"""

    rlm7w_panel = ['15473 Zygosity Call', '62768 Zygosity Call', '295723 Zygosity Call']

    #  index the headers to see if the CR panel is in the kraken file, and if so, which column #'s they are in
    rlm7w_index = []
    for rlm7wi in rlm7w_panel:
        # noinspection PyBroadException
        try:
            rlm7w_index.append(kheaders.index(rlm7wi))
        except:
            continue
    rlm7w_len = int(len(rlm7w_index))

    # if that index is not empty, make a pandas dataframe, and write that section to a new csv to run analysis
    if rlm7w_len >= 2:
        # keep the unique identifiers from kraken study in new temp dataframe
        rlm7w_headers = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']

        for i in rlm7w_index:  # add the RLM7 winter assays in the kraken file to the header list
            rlm7w_headers.append(kheaders[i])

        rlm7wdf = kdf[rlm7w_headers]  # make a dataframe of just RLM7 winter data
        rlm7wdf.to_csv('rlm7w_temp.csv', index=False)

        with open('rlm7w_temp.csv', 'rb') as rlm7winput:  # open that newly created dataframe
            rlm7wreader = csv.DictReader(rlm7winput)
            rlm7wheaders = rlm7wreader.fieldnames

            with open('rlm7w_call.csv', 'wb') as rlm7woutput:  # make a new file with the rlm7 call column added
                rlm7wheaders.append('RLM 7 Winter Summary')
                rlm7wwriter = csv.DictWriter(rlm7woutput, fieldnames=rlm7wheaders)
                rlm7wwriter.writeheader()

                for line in rlm7wreader:

                    dict_rlm7w = {}  # make a dictionary for each line in temp file, to pass into the TLM class
                    if len(rlm7w_index) > 0:
                        for i in rlm7w_index:
                            dict_rlm7w[kheaders[i]] = line[kheaders[i]]

                    rlm7w = Trait(dict_rlm7w)  # initiate the TLM class with dictionary
                    rlm7w_call = rlm7w.get_tlm_call()  # get the trait/seg/wildtype call for this line in the file
                    line['RLM 7 Winter Summary'] = rlm7w_call  # add call to the DictReader line

                    rlm7wwriter.writerow(line)  # write the line (with the call now in the line) to the temp call file

                pandas_list.append('rlm7w')  # add rlm7w (string) to pandas list for prompts to merge later

        os.remove('rlm7w_temp.csv')

    # ---------------------------------------------- BL10 Analysis ----------------------------------------------------
    """determine if BL10 panel is in the kraken study, and if so makes a df with a BL10 call"""

    bl10_panel = ['302267 Zygosity Call', '331689 Zygosity Call', '281051 Zygosity Call', '303572 Zygosity Call',
                  '303605 Zygosity Call', '303806 Zygosity Call']

    #  index the headers to see if the CR panel is in the kraken file, and if so, which column #'s they are in
    bl10_index = []
    for bl10i in bl10_panel:
        # noinspection PyBroadException
        try:
            bl10_index.append(kheaders.index(bl10i))
        except:
            continue
    bl10_len = int(len(bl10_index))

    # if that index is not empty, make a pandas dataframe, and write that section to a new csv to run analysis
    if bl10_len >= 2:
        #  keep the unique identifiers from kraken study in new temp dataframe
        bl10_headers = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']

        for i in bl10_index:  # add the BL 10 assays in the kraken file to the header list
            bl10_headers.append(kheaders[i])

        bl10df = kdf[bl10_headers]  # make a dataframe of just BL 10 data
        bl10df.to_csv('bl10_temp.csv', index=False)

        with open('bl10_temp.csv', 'rb') as bl10input:  # open that newly created dataframe
            bl10reader = csv.DictReader(bl10input)
            bl10headers = bl10reader.fieldnames

            with open('bl10_call.csv', 'wb') as bl10output:  # make a new file with the CR call column added
                bl10headers.append('BL10 Summary')
                bl10writer = csv.DictWriter(bl10output, fieldnames=bl10headers)
                bl10writer.writeheader()

                for line in bl10reader:

                    dict_bl10 = {}  # make a dictionary for each line in temp file, to pass into the TLM class
                    if len(bl10_index) > 0:
                        for i in bl10_index:
                            dict_bl10[kheaders[i]] = line[kheaders[i]]

                    bl10 = Trait(dict_bl10)  # initiate the TLM class with the dictionary
                    bl10_call = bl10.get_tlm_call()  # get the trait/seg/wildtype call for this line in the file
                    line['BL10 Summary'] = bl10_call  # add call to the DictReader line

                    bl10writer.writerow(line)  # write the line (with the call now in the line) to the temp call file

                pandas_list.append('bl10')  # add bl10 (string) to pandas list for prompts to merge later

        os.remove('bl10_temp.csv')

    # ---------------------------------------------- N13 Analysis ----------------------------------------------------
    """determine if N13 panel is in the kraken study, and if so makes a df with a N13 call"""

    n13_panel = ['15925 Zygosity Call', '106265 Zygosity Call', '220074 Zygosity Call', '320634 Zygosity Call',
                 '320752 Zygosity Call', '320826 Zygosity Call', '328092 Zygosity Call', '332445 Zygosity Call',
                 '335227 Zygosity Call', 'N13-7444161 Zygosity Call', 'N13-7444297 Zygosity Call',
                 'N13-7444437 Zygosity Call']

    #  index the headers to see if the N13 panel is in the kraken file, and if so, which column #'s they are in
    n13_index = []
    for n13i in n13_panel:
        # noinspection PyBroadException
        try:
            n13_index.append(kheaders.index(n13i))
        except:
            continue
    n13_len = int(len(n13_index))

    # if that index is not empty, make a pandas dataframe, and write that section to a new csv to run analysis
    if n13_len >= 2:
        #  keep the unique identifiers from kraken study in new temp dataframe
        n13_headers = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']

        for i in n13_index:  # add the N13 assays in the kraken file to the header list
            n13_headers.append(kheaders[i])

        n13df = kdf[n13_headers]  # make a dataframe of just N13 data
        n13df.to_csv('n13_temp.csv', index=False)

        with open('n13_temp.csv', 'rb') as n13input:  # open that newly created dataframe
            n13reader = csv.DictReader(n13input)
            n13headers = n13reader.fieldnames

            with open('n13_call.csv', 'wb') as n13output:  # make a new file with the CR call column added
                n13headers.append('ACM N13 Summary')
                n13writer = csv.DictWriter(n13output, fieldnames=n13headers)
                n13writer.writeheader()

                for line in n13reader:

                    dict_n13 = {}  # make a dictionary for each line in temp file, to pass into the TLM class
                    if len(n13_index) > 0:
                        for i in n13_index:
                            dict_n13[kheaders[i]] = line[kheaders[i]]

                    n13 = Trait(dict_n13)  # initiate the TLM class with the dictionary
                    n13_call = n13.get_tlm_call()  # get the trait/seg/wildtype call for this line in the file
                    line['ACM N13 Summary'] = n13_call  # add call to the DictReader line

                    n13writer.writerow(line)  # write the line (with the call now in the line) to the temp call file

                pandas_list.append('n13')  # add n13 (string) to pandas list for prompts to merge later

        os.remove('n13_temp.csv')

        # -------------------------------------------- FAE Analysis -------------------------------------------------
        """determine if FAE panel is in the kraken study, and if so makes a df with a FAE summary call"""

        fae_panel = ['FAE 1-1 A Zygosity Call', 'FAE 1-2 C Zygosity Call']

        #  index the headers to see if the FAE panel is in the kraken file, and if so, which column #'s they are in
        fae_index = []
        for faei in fae_panel:
            # noinspection PyBroadException
            try:
                fae_index.append(kheaders.index(faei))
            except:
                continue
        fae_len = int(len(fae_index))

        # if that index is not empty, make a pandas dataframe, and write that section to a new csv to run analysis
        if fae_len >= 2:
            #  keep the unique identifiers from kraken study in new temp dataframe
            fae_headers = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']

            for i in fae_index:  # add the FAE assays in the kraken file to the header list
                fae_headers.append(kheaders[i])

            fae_df = kdf[fae_headers]  # make a dataframe of just FAE data
            fae_df.to_csv('fae_temp.csv', index=False)

            with open('fae_temp.csv', 'rb') as faeinput:  # open that newly created dataframe
                faereader = csv.DictReader(faeinput)
                faeheaders = faereader.fieldnames

                with open('fae_call.csv', 'wb') as faeoutput:  # make a new file with the FAE call column added
                    faeheaders.append('FAE Summary Call')
                    faewriter = csv.DictWriter(faeoutput, fieldnames=faeheaders)
                    faewriter.writeheader()

                    for line in faereader:

                        dict_fae = {}  # make a dictionary for each line in temp file, to pass into the TLM class
                        if len(fae_index) > 0:
                            for i in fae_index:
                                dict_fae[kheaders[i]] = line[kheaders[i]]

                        fae = Trait(dict_fae)  # initiate the TLM class with the dictionary
                        fae_call = fae.get_fae()  # get the trait/seg/wildtype call for this line in the file
                        line['FAE Summary Call'] = fae_call  # add call to the DictReader line

                        faewriter.writerow(line)  # write the line (with the call now in the line) to the call file

                    pandas_list.append('fae')  # add n13 (string) to pandas list for prompts to merge later

            os.remove('fae_temp.csv')

    # ----------------------------------------- Merge with Pandas -----------------------------------------------------

    # read each panels analyzed data as a dataframe, and merge them in sequence.

    krak_info = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']

    df_initial = df_final[krak_info]  # sample and well info dataframe only (no molecular data)

    # merge initial with N9 dataframe if panel is in the kraken file

    if 'n9' in pandas_list:
        n9_df = pd.read_csv('n9_call.csv')
        merge1 = pd.merge(left=df_initial, right=n9_df, how='left')
        os.remove('n9_call.csv')
    else:
        merge1 = df_initial

    # merge that with n13 dataframe if panel is in the kraken file

    if 'n13' in pandas_list:
        n13_df = pd.read_csv('n13_call.csv')
        merge2 = pd.merge(left=merge1, right=n13_df, how='left')
        os.remove('n13_call.csv')
    else:
        merge2 = merge1

    # merge that with castle dataframe if panel is in the kraken file

    if 'castle' in pandas_list:
        cas_df = pd.read_csv('cas_call.csv')
        merge3 = pd.merge(left=merge2, right=cas_df, how='left')
        os.remove('cas_call.csv')
    else:
        merge3 = merge2

    # merge that with clubroot dataframe if panel is in the kraken file

    if 'cr' in pandas_list:
        cr_df = pd.read_csv('cr_call.csv')
        merge4 = pd.merge(left=merge3, right=cr_df, how='left')
        os.remove('cr_call.csv')
    else:
        merge4 = merge3

    # merge that with RLM7 spring dataframe if panel is in the kraken file

    if 'rlm7s' in pandas_list:
        rlm7s_df = pd.read_csv('rlm7s_call.csv')
        merge5 = pd.merge(left=merge4, right=rlm7s_df, how='left')
        os.remove('rlm7s_call.csv')
    else:
        merge5 = merge4

    # merge that with RLM7 winter dataframe if that panel is in the kraken file

    if 'rlm7w' in pandas_list:
        rlm7w_df = pd.read_csv('rlm7w_call.csv')
        merge6 = pd.merge(left=merge5, right=rlm7w_df, how='left')
        os.remove('rlm7w_call.csv')
    else:
        merge6 = merge5

    # merge that with BL10 dataframe if that panel is in the kraken file

    if 'bl10' in pandas_list:
        bl10_df = pd.read_csv('bl10_call.csv')
        merge7 = pd.merge(left=merge6, right=bl10_df, how='left')
        os.remove('bl10_call.csv')
    else:
        merge7 = merge6

    # merge that with FAE dataframe if that panel is in the kraken file

    if 'fae' in pandas_list:
        fae_df = pd.read_csv('fae_call.csv')
        merge8 = pd.merge(left=merge7, right=fae_df, how='left')
        os.remove('fae_call.csv')
    else:
        merge8 = merge7

    # merge with final dataframe and write the file (gets any assays not in panel tacked on to end of kraken report)

    merge_final = pd.merge(left=merge8, right=df_final, how='left')

    os.chdir(completed_path)
    merge_final.to_excel(results_writer, 'Report', index=False)
    summary_df.to_excel(results_writer, 'Summary Table')
    results_writer.save()

    os.chdir(path)  # change back to root folder
    os.remove(encodename)

    print 'Successfully completed ' + fname + '...\n'

    time.sleep(1)

    os.remove(fname)

print 'Analysis completed on all files'

time.sleep(2)
