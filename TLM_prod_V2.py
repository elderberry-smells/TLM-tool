import csv
import pandas as pd
from pandas import ExcelWriter
from mas_class_tools import *
from summary_class import *
import glob
import os
import time

'''Version 2 of the TLM panel tool-- Validated March 2017 --Author u590135'''

#  get a list of all the csv files to run the script on
path = r'C:\Users\u590135\Code\testing code\Trait and Seg Tool'
extension = 'csv'
os.chdir(path)
result = [i for i in glob.glob('*.{}'.format(extension))]

completed_path = r'C:\Users\u590135\Code\testing code\Trait and Seg Tool\completed'

results_files = []
for i in result:
    if 'temp_' in i:
        continue
    if '_call' in i:
        continue
    else:
        results_files.append(i)


# ---------------------------------------- Defined Variables ----------------------------------------------------------


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


def cyto_gt200(header_list):
    """if header list from the Kraken file have either A Cyto or GT200 change the listed name to the Variety name
    :param header_list: is a list of headers found in the Kraken file being analyzed """

    if 'Cyto Zygosity Call' in header_list:
        #  change the header to read Cyto Call
        header_list = ['Cyto Call' if hdr == 'Cyto Zygosity Call' else hdr for hdr in header_list]

    if 'GT200 Zygosity Call' in header_list:
        #  change the header to read GT200 Call
        header_list = ['GT200 Call' if hdr2 == 'GT200 Zygosity Call' else hdr2 for hdr2 in header_list]

    return header_list


def cgl_conv(header_list):
    """get a list of what headers are in the kraken study, and assign values to cyto, gt200 and LepR3
    :param header_list: header list from Kraken csv file
    """

    cgl_sum = 0
    if 'Cyto Zygosity Call' in header_list:
        cgl_sum += 1

    if 'GT 200 Zygosity Call' in header_list:
        cgl_sum += 10

    if 'LepR3-R3S3A Zygosity Call' in header_list:
        cgl_sum += 100

    return cgl_sum


#  ----------------------------------  read the folder and find files for analysis  -----------------------------------

for fname in results_files:
    print('Analyzing ' + fname + '...')

    #  make the naming conventions for the files being written in the end
    encodename = 'encoded_' + fname
    find_the_dot = fname.find('.')
    final_name = fname[:find_the_dot] + '_completed.xlsx'
    results_writer = ExcelWriter(final_name, options={'encoding': 'utf-8'})

    #  ---------------------------------------  Encoding the file to UTF-8  -------------------------------------------
    with open(fname, 'r+') as encode_kraken:
        encode_reader = csv.DictReader(encode_kraken)
        encode_headers = encode_reader.fieldnames

    with open(encodename, 'wb') as encode_file:  # write a temporary UTF-8 encoded version of the kraken file
        writer = csv.DictWriter(encode_file, fieldnames=encode_headers)

        r = open(fname, 'rb')

        for line in r:
            bits = line.split(',')
            bits[3] = bits[3].decode('latin-1').encode('utf-8')  # encode the pedigree line as UTF-8 instead of Latin-1
            encode_file.write(','.join(bits))

        r.close()

    #  ----------------------------------------  Summary Table  --------------------------------------------------------

    with open(encodename, 'rb') as summary_kraken:
        sum_reader = csv.DictReader(summary_kraken)
        sum_headers = sum_reader.fieldnames
        sum_assays = [i for i in sum_headers if 'Zygosity' in i]

        summary_init = SummaryTable(sum_reader, sum_headers)
        summary_dict = summary_init.get_summary()
        summary_df1 = pd.DataFrame.from_dict(summary_dict)
        summary_df2 = summary_df1[sum_assays]
        summary_df = summary_df2.reindex(["Trait", "Seg", "Wildtype", "No Call", "Fail", "% Data Return"])

    # ------------------------------------------  Analysis  ----------------------------------------------------------

    #  Determine what assays were run in the Kraken study, and if anything needs to be analyzed or converted in the file
    with open(encodename, 'rb') as kraken:
        kreader = csv.DictReader(kraken)
        kheaders = kreader.fieldnames
        new_headers = cyto_gt200(kheaders)
        conv_kraken = cgl_conv(kheaders)
        creader = csv.reader(kraken)

        #  ----------------------------------- Cyto, GT200, and LepR3 Conversion --------------------------------------

        if conv_kraken > 0:  # if any of the three assays are in the kraken study.

            # get the index of which columns they are in
            cyto = ['Cyto Zygosity Call']
            gt200 = ['GT 200 Zygosity Call']
            lep3 = ['LepR3-R3S3A Zygosity Call']

            cyto_index = []
            for ci in cyto:
                # noinspection PyBroadException
                try:
                    cyto_index.append(kheaders.index(ci))
                except:
                    continue

            gt200_index = []
            for gi in gt200:
                # noinspection PyBroadException
                try:
                    gt200_index.append(kheaders.index(gi))
                except:
                    continue

            lep3_index = []
            for li in lep3:
                # noinspection PyBroadException
                try:
                    lep3_index.append(kheaders.index(li))
                except:
                    continue

            # re-write headers so that if Acyto or GT200 are in list, they are overwritten by new column names
            with open('temp_' + fname, 'wb') as outfile:
                hdrwriter = csv.writer(outfile)
                hdrwriter.writerow(new_headers)

                for line in creader:

                    if conv_kraken == 1:  # Kraken file has A Cyto only
                        cyto_col = int(cyto_index[0])  # get the A cyto column number from the index
                        cy = CallConversion(line[cyto_col])  # call on class CallConverison to convert cyto row data
                        cy_call = cy.cyto_call()  # return the call as A-Cyto or B-Cyto
                        line[cyto_col] = line[cyto_col].replace(line[cyto_col], cy_call)  # replace call
                        hdrwriter.writerow(line)  # write the line to the new temporary file

                    elif conv_kraken == 10:  # Kraken file has GT200 only
                        gt_col = int(gt200_index[0])
                        gt = CallConversion(line[gt_col])
                        gt_call = gt.plus_minus()
                        line[gt_col] = line[gt_col].replace(line[gt_col], gt_call)
                        hdrwriter.writerow(line)

                    elif conv_kraken == 11:  # Kraken file has ACyto and GT200
                        cyto_col = int(cyto_index[0])
                        gt_col = int(gt200_index[0])
                        cy = CallConversion(line[cyto_col])
                        cy_call = cy.cyto_call()
                        gt = CallConversion(line[gt_col])
                        gt_call = gt.plus_minus()
                        line[cyto_col] = line[cyto_col].replace(line[cyto_col], cy_call)
                        line[gt_col] = line[gt_col].replace(line[gt_col], gt_call)
                        hdrwriter.writerow(line)

                    elif conv_kraken == 100:  # Kraken file has LepR3 only
                        l3_col = int(lep3_index[0])
                        l3 = CallConversion(line[l3_col])
                        l3_call = l3.zygo_conv()
                        line[l3_col] = line[l3_col].replace(line[l3_col], l3_call)
                        hdrwriter.writerow(line)

                    elif conv_kraken == 101:  # Kraken file has A Cyto and LepR3
                        cyto_col = int(cyto_index[0])
                        l3_col = int(lep3_index[0])
                        cy = CallConversion(line[cyto_col])
                        cy_call = cy.cyto_call()
                        l3 = CallConversion(line[l3_col])
                        l3_call = l3.zygo_conv()
                        line[cyto_col] = line[cyto_col].replace(line[cyto_col], cy_call)
                        line[l3_col] = line[l3_col].replace(line[l3_col], l3_call)
                        hdrwriter.writerow(line)

                    elif conv_kraken == 110:  # Kraken file has GT200 and LepR3
                        gt_col = int(gt200_index[0])
                        l3_col = int(lep3_index[0])
                        gt = CallConversion(line[gt_col])
                        gt_call = gt.plus_minus()
                        l3 = CallConversion(line[l3_col])
                        l3_call = l3.zygo_conv()
                        line[gt_col] = line[gt_col].replace(line[gt_col], gt_call)
                        line[l3_col] = line[l3_col].replace(line[l3_col], l3_call)
                        hdrwriter.writerow(line)

                    elif conv_kraken == 111:  # Kraken file has A Cyto, GT200, and LepR3
                        cyto_col = int(cyto_index[0])
                        gt_col = int(gt200_index[0])
                        l3_col = int(lep3_index[0])
                        cy = CallConversion(line[cyto_col])
                        cy_call = cy.cyto_call()
                        gt = CallConversion(line[gt_col])
                        gt_call = gt.plus_minus()
                        l3 = CallConversion(line[l3_col])
                        l3_call = l3.zygo_conv()
                        line[cyto_col] = line[cyto_col].replace(line[cyto_col], cy_call)
                        line[gt_col] = line[gt_col].replace(line[gt_col], gt_call)
                        line[l3_col] = line[l3_col].replace(line[l3_col], l3_call)
                        hdrwriter.writerow(line)

            kdf = pd.read_csv('temp_' + fname)
            df_final = pd.read_csv('temp_' + fname)
            os.remove('temp_' + fname)

        else:
            kdf = pd.read_csv(encodename)
            df_final = pd.read_csv(encodename)

    pandas_list = []  # an empty list to append each panels dataframe name for a final merge sequence

    # ----------------------------------------- N9 Panel --------------------------------------------------------------

    """determine if N9 panel is in the kraken study, and if so makes a df with a N9 call"""

    n9_panel = ['DBSNP357223 Zygosity Call', 'DBSNP357253 Zygosity Call', 'DBSNP357255 Zygosity Call',
                'DBSNP357256 Zygosity Call']

    #  index the headers to see if the N9 panel is in the kraken file, and if so, which column #'s they are in
    n9_index = []
    for i in n9_panel:
        # noinspection PyBroadException
        try:
            n9_index.append(kheaders.index(i))
        except:
            continue
    n9_len = int(len(n9_index))

    # if N9 panel is in the kraken file, make a pandas dataframe, and write that section to a new csv to run analysis
    if n9_len == 4:
        n9_headers = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']
        for i in n9_index:
            n9_headers.append(kheaders[i])
        n9df = kdf[n9_headers]
        n9df.to_csv('n9_temp.csv', index=False)

        with open('n9_temp.csv', 'rb') as n9input:
            n9reader = csv.DictReader(n9input)
            n9headers = n9reader.fieldnames
            n9headers.append('ACM N9 Summary')

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
                                       'ACM N9 Summary': n9_call})

                pandas_list.append('n9')

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
    if cas_len == 4:
        castle_headers = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']
        for i in castle_index:
            castle_headers.append(kheaders[i])
        casdf = kdf[castle_headers]
        casdf.to_csv('cas_temp.csv', index=False)

        with open('cas_temp.csv', 'rb') as casinput:
            csreader = csv.DictReader(casinput)
            csheaders = csreader.fieldnames
            csheaders.append('Castle Summary')

            with open('cas_call.csv', 'wb') as casoutput:
                cswriter = csv.DictWriter(casoutput, fieldnames=csheaders)
                cswriter.writeheader()

                for line in csreader:
                    casid = TraitCastle(line['298512 Zygosity Call'], line['298516 Zygosity Call'],
                                        line['298518 Zygosity Call'], line['298535 Zygosity Call'])
                    cas_call = casid.cas_call()

                    cswriter.writerow({'Box': line['Box'], 'Well': line['Well'], 'Project': line['Project'],
                                       'Pedigree': line['Pedigree'], 'Source ID': line['Source ID'],
                                       'Geno_Id': line['Geno_Id'], 'RowId': line['RowId'], 'Loc Seq#': line['Loc Seq#'],
                                       '298512 Zygosity Call': line['298512 Zygosity Call'],
                                       '298516 Zygosity Call': line['298516 Zygosity Call'],
                                       '298518 Zygosity Call': line['298518 Zygosity Call'],
                                       '298535 Zygosity Call': line['298535 Zygosity Call'],
                                       'Castle Summary': cas_call})

                pandas_list.append('castle')

        os.remove('cas_temp.csv')

    # ----------------------------------------- CR Panel --------------------------------------------------------------
    """determine if Castle panel is in the kraken study, and if so makes a df with a Castle call"""

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
    if cr_len == 2:
        cr_headers = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']
        for i in cr_index:
            cr_headers.append(kheaders[i])
        crdf = kdf[cr_headers]
        crdf.to_csv('cr_temp.csv', index=False)

        with open('cr_temp.csv', 'rb') as crinput:
            crreader = csv.DictReader(crinput)
            crheaders = crreader.fieldnames
            crheaders.append('CR Mendel Summary')

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
                                       'CR Mendel Summary': cr_call})

                pandas_list.append('cr')

        os.remove('cr_temp.csv')
    # ----------------------------------------- Merge with Pandas -----------------------------------------------------

    # read each panels analyzed data as a dataframe, and merge them in sequence.

    krak_info = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']

    df_initial = df_final[krak_info]  # sample and well info (no molecular data)

    merge_number = determine_merge(pandas_list)

    #  begin merging in sequence. start with initial df, merging the available df's, and end with merge on the final df

    if merge_number == 1:  # if the panel consists of only N9
        n9_df = pd.read_csv('n9_call.csv')

        merge_init = pd.merge(left=df_initial, right=n9_df, how='left')
        merge_final = pd.merge(left=merge_init, right=df_final, how='left')
        os.remove('n9_call.csv')

        os.chdir(completed_path)  # change to path for completed files
        merge_final.to_excel(results_writer, 'Report', index=False)
        summary_df.to_excel(results_writer, 'Summary Table')
        results_writer.save()
        os.chdir(path)  # change back to root folder
        os.remove(encodename)

    elif merge_number == 10:  # if the panel consists of only Castle
        cas_df = pd.read_csv('cas_call.csv')

        merge_init = pd.merge(left=df_initial, right=cas_df, how='left')
        merge_final = pd.merge(left=merge_init, right=df_final, how='left')
        os.remove('cas_call.csv')
        os.chdir(completed_path)  # change to path for completed files
        merge_final.to_excel(results_writer, 'Report', index=False)
        summary_df.to_excel(results_writer, 'Summary Table')
        results_writer.save()
        os.chdir(path)  # change back to root folder
        os.remove(encodename)

    elif merge_number == 11:  # if the panel consists of N9 and Castle
        n9_df = pd.read_csv('n9_call.csv')
        cas_df = pd.read_csv('cas_call.csv')

        merge_init = pd.merge(left=df_initial, right=n9_df, how='left')
        merge_two = pd.merge(left=merge_init, right=cas_df, how='left')
        merge_final = pd.merge(left=merge_two, right=df_final, how='left')
        os.remove('n9_call.csv')
        os.remove('cas_call.csv')
        os.chdir(completed_path)  # change to path for completed files
        merge_final.to_excel(results_writer, 'Report', index=False)
        summary_df.to_excel(results_writer, 'Summary Table')
        results_writer.save()
        os.chdir(path)  # change back to root folder
        os.remove(encodename)

    elif merge_number == 100:  # if the panel consists of only CR
        cr_df = pd.read_csv('cr_call.csv')

        merge_init = pd.merge(left=df_initial, right=cr_df, how='left')
        merge_final = pd.merge(left=merge_init, right=df_final, how='left')
        os.remove('cr_call.csv')
        os.chdir(completed_path)  # change to path for completed files
        merge_final.to_excel(results_writer, 'Report', index=False)
        summary_df.to_excel(results_writer, 'Summary Table')
        results_writer.save()
        os.chdir(path)  # change back to root folder
        os.remove(encodename)

    elif merge_number == 101:  # if the panel consists of N9 and CR
        n9_df = pd.read_csv('n9_call.csv')
        cr_df = pd.read_csv('cr_call.csv')

        merge_init = pd.merge(left=df_initial, right=n9_df, how='left')
        merge_two = pd.merge(left=merge_init, right=cr_df, how='left')
        merge_final = pd.merge(left=merge_two, right=df_final, how='left')
        os.remove('n9_call.csv')
        os.remove('cr_call.csv')
        os.chdir(completed_path)  # change to path for completed files
        merge_final.to_excel(results_writer, 'Report', index=False)
        summary_df.to_excel(results_writer, 'Summary Table')
        results_writer.save()
        os.chdir(path)  # change back to root folder
        os.remove(encodename)

    elif merge_number == 110:  # if the panel consists of Castle and CR
        cas_df = pd.read_csv('cas_call.csv')
        cr_df = pd.read_csv('cr_call.csv')

        merge_init = pd.merge(left=df_initial, right=cas_df, how='left')
        merge_two = pd.merge(left=merge_init, right=cr_df, how='left')
        merge_final = pd.merge(left=merge_two, right=df_final, how='left')
        os.remove('cas_call.csv')
        os.remove('cr_call.csv')
        os.chdir(completed_path)  # change to path for completed files
        merge_final.to_excel(results_writer, 'Report', index=False)
        summary_df.to_excel(results_writer, 'Summary Table')
        results_writer.save()
        os.chdir(path)  # change back to root folder
        os.remove(encodename)

    elif merge_number == 111:  # if the panel consists of N9, Castle, and CR
        n9_df = pd.read_csv('n9_call.csv')
        cas_df = pd.read_csv('cas_call.csv')
        cr_df = pd.read_csv('cr_call.csv')

        merge_init = pd.merge(left=df_initial, right=n9_df, how='left')
        merge_two = pd.merge(left=merge_init, right=cas_df, how='left')
        merge_three = pd.merge(left=merge_two, right=cr_df, how='left')
        merge_final = pd.merge(left=merge_three, right=df_final, how='left')
        os.remove('n9_call.csv')
        os.remove('cas_call.csv')
        os.remove('cr_call.csv')
        os.chdir(completed_path)  # change to path for completed files
        merge_final.to_excel(results_writer, 'Report', index=False)
        summary_df.to_excel(results_writer, 'Summary Table')
        results_writer.save()
        os.chdir(path)  # change back to root folder
        os.remove(encodename)

    else:
        os.chdir(completed_path)  # change to path for completed files
        df_final.to_excel(results_writer, 'Report', index=False)
        summary_df.to_excel(results_writer, 'Summary Table')
        results_writer.save()
        os.chdir(path)  # change back to root folder
        os.remove(encodename)

    print 'Successfully completed ' + fname + '...\n'

    time.sleep(2)

print 'Analysis completed on all files'

time.sleep(2)
