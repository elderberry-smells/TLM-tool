class SummaryTable(object):
    """grab the reader data from kraken file and make a dictionary of the # of calls for each assay"""

    def __init__(self, kraken, headers):
        self.kraken = kraken
        self.headers = [i for i in headers if 'Zygosity' in i]  # only grab headers that have 'zygosity' in name

    def get_summary(self):
        """update counts of each assays calls and return the dictionary"""

        assay_dic = {}  # empty dictionary to append assays and calls to

        # Trait = G:G, Seg = G:A, Wildtype = A:A
        type_1 = [
            'DBSNP357223 Zygosity Call',
            'DBSNP357253 Zygosity Call',
            'DBSNP357256 Zygosity Call',
            '298535 Zygosity Call',
            '295699 Zygosity Call',
            '303605 Zygosity Call',
            'DBSNP357252 Zygosity Call',
            'N13-7444437 Zygosity Call'
        ]

        # Trait = G:G, Seg = G:T, Wildtype = T:T
        type_2 = [
            '298512 Zygosity Call'
        ]

        # Trait = T:T, Seg = G:T, Wildtype = G:G
        type_3 = [
            '298516 Zygosity Call',
            '303806 Zygosity Call'
        ]

        # Trait = T:T, Seg = C:T, Wildtype = C:C
        type_4 = [
            '298518 Zygosity Call',
            'DBSNP357255 Zygosity Call',
            '295698 Zygosity Call',
            '302267 Zygosity Call',
            '281051 Zygosity Call',
            '220074 Zygosity Call',
            '320752 Zygosity Call'
        ]

        # Trait = T:T, Seg = T:A, Wildtype = A:A
        type_5 = [
            '15473 Zygosity Call'
        ]

        # Trait = A:A, Seg = G:A, Wildtype = G:G
        type_6 = [
            '62768 Zygosity Call',
            '295723 Zygosity Call',
            'N13-7444161 Zygosity Call',
            '335227 Zygosity Call'
        ]

        # Trait = C:C, Seg = C:A, Wildtype = A:A
        type_7 = [
            '331689 Zygosity Call',
            '303572 Zygosity Call',
            'N13-7444297 Zygosity Call'
        ]

        # Trait = G:G, Seg = G:C, Wildtype = C:C
        type_8 = [
            'DBSNP357226 Zygosity Call',
            '15925 Zygosity Call'
        ]

        # Trait = C:C, Seg = C:T, Wildtype = T:T
        type_9 = [
            'DBSNP357233 Zygosity Call',
            'DBSNP357257 Zygosity Call',
            '332445 Zygosity Call',
            '320634 Zygosity Call'
        ]

        # Trait = A:A, Seg = T:A, Wildtype = 'T:T
        type_10 = [
            '320826 Zygosity Call',
            '106265 Zygosity Call'
        ]

        # Trait = A:A, Seg = C:A, Wildtype = C:C
        type_11 = [
            '328092 Zygosity Call'
        ]

        #  append all assays in kraken file to the dictionary with 0 as blank call to begin with
        for i in self.headers:
            assay_dic[i] = {'Trait': 0, 'Seg': 0, 'Wildtype': 0, 'No Call': 0, 'Fail': 0}

        #  read the kraken study line by line and start counting the homo/hemi/null values in the dictionary for
        #  each assay
        for row in self.kraken:

            for assay in self.headers:

                if assay in type_1:  # if a SNP call from specific set 1, set homo/hemi/null to that grouping
                    if row[assay] == 'G:G':
                        assay_dic[assay]['Trait'] += 1
                    elif row[assay] == 'G:A':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'A:G':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'A:A':
                        assay_dic[assay]['Wildtype'] += 1
                    elif row[assay] == 'No Call':
                        assay_dic[assay]['No Call'] += 1
                    elif row[assay] == 'No Data':
                        assay_dic[assay]['Fail'] += 1
                    else:
                        continue

                elif assay in type_2:  # if a SNP call from specific set 2, set homo/hemi/null to that grouping
                    if row[assay] == 'G:G':
                        assay_dic[assay]['Trait'] += 1
                    elif row[assay] == 'G:T':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'T:G':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'T:T':
                        assay_dic[assay]['Wildtype'] += 1
                    elif row[assay] == 'No Call':
                        assay_dic[assay]['No Call'] += 1
                    elif row[assay] == 'No Data':
                        assay_dic[assay]['Fail'] += 1
                    else:
                        continue

                elif assay in type_3:  # if a SNP call from specific set 3, set homo/hemi/null to that grouping
                    if row[assay] == 'T:T':
                        assay_dic[assay]['Trait'] += 1
                    elif row[assay] == 'G:T':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'T:G':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'G:G':
                        assay_dic[assay]['Wildtype'] += 1
                    elif row[assay] == 'No Call':
                        assay_dic[assay]['No Call'] += 1
                    elif row[assay] == 'No Data':
                        assay_dic[assay]['Fail'] += 1
                    else:
                        continue

                elif assay in type_4:  # if a SNP call from specific set 4, set homo/hemi/null to that grouping
                    if row[assay] == 'T:T':
                        assay_dic[assay]['Trait'] += 1
                    elif row[assay] == 'C:T':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'T:C':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'C:C':
                        assay_dic[assay]['Wildtype'] += 1
                    elif row[assay] == 'No Call':
                        assay_dic[assay]['No Call'] += 1
                    elif row[assay] == 'No Data':
                        assay_dic[assay]['Fail'] += 1
                    else:
                        continue

                elif assay in type_5:  # if a SNP call from specific set 5, set homo/hemi/null to that grouping
                    if row[assay] == 'T:T':
                        assay_dic[assay]['Trait'] += 1
                    elif row[assay] == 'T:A':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'A:T':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'A:A':
                        assay_dic[assay]['Wildtype'] += 1
                    elif row[assay] == 'No Call':
                        assay_dic[assay]['No Call'] += 1
                    elif row[assay] == 'No Data':
                        assay_dic[assay]['Fail'] += 1
                    else:
                        continue

                elif assay in type_6:  # if a SNP call from specific set 6, set homo/hemi/null to that grouping
                    if row[assay] == 'A:A':
                        assay_dic[assay]['Trait'] += 1
                    elif row[assay] == 'G:A':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'A:G':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'G:G':
                        assay_dic[assay]['Wildtype'] += 1
                    elif row[assay] == 'No Call':
                        assay_dic[assay]['No Call'] += 1
                    elif row[assay] == 'No Data':
                        assay_dic[assay]['Fail'] += 1
                    else:
                        continue

                elif assay in type_7:  # if a SNP call from specific set 7, set homo/hemi/null to that grouping
                    if row[assay] == 'C:C':
                        assay_dic[assay]['Trait'] += 1
                    elif row[assay] == 'C:A':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'A:C':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'A:A':
                        assay_dic[assay]['Wildtype'] += 1
                    elif row[assay] == 'No Call':
                        assay_dic[assay]['No Call'] += 1
                    elif row[assay] == 'No Data':
                        assay_dic[assay]['Fail'] += 1
                    else:
                        continue

                elif assay in type_8:  # if a SNP call from specific set 8, set homo/hemi/null to that grouping
                    if row[assay] == 'G:G':
                        assay_dic[assay]['Trait'] += 1
                    elif row[assay] == 'G:C':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'C:G':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'C:C':
                        assay_dic[assay]['Wildtype'] += 1
                    elif row[assay] == 'No Call':
                        assay_dic[assay]['No Call'] += 1
                    elif row[assay] == 'No Data':
                        assay_dic[assay]['Fail'] += 1
                    else:
                        continue

                elif assay in type_9:  # if a SNP call from specific set 9, set homo/hemi/null to that grouping
                    if row[assay] == 'C:C':
                        assay_dic[assay]['Trait'] += 1
                    elif row[assay] == 'C:T':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'T:C':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'T:T':
                        assay_dic[assay]['Wildtype'] += 1
                    elif row[assay] == 'No Call':
                        assay_dic[assay]['No Call'] += 1
                    elif row[assay] == 'No Data':
                        assay_dic[assay]['Fail'] += 1
                    else:
                        continue

                elif assay in type_10:  # if a SNP call from specific set 10, set homo/hemi/null to that grouping
                    if row[assay] == 'A:A':
                        assay_dic[assay]['Trait'] += 1
                    elif row[assay] == 'T:A':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'A:T':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'T:T':
                        assay_dic[assay]['Wildtype'] += 1
                    elif row[assay] == 'No Call':
                        assay_dic[assay]['No Call'] += 1
                    elif row[assay] == 'No Data':
                        assay_dic[assay]['Fail'] += 1
                    else:
                        continue

                elif assay in type_11:  # if a SNP call from specific set 11, set homo/hemi/null to that grouping
                    if row[assay] == 'A:A':
                        assay_dic[assay]['Trait'] += 1
                    elif row[assay] == 'C:A':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'A:C':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'C:C':
                        assay_dic[assay]['Wildtype'] += 1
                    elif row[assay] == 'No Call':
                        assay_dic[assay]['No Call'] += 1
                    elif row[assay] == 'No Data':
                        assay_dic[assay]['Fail'] += 1
                    else:
                        continue

                else:  # if assays are homo/hemi/null or trait/seg/null run the counter
                    if row[assay] == 'Homo':
                        assay_dic[assay]['Trait'] += 1
                    elif row[assay] == 'Trait':
                        assay_dic[assay]['Trait'] += 1

                    elif row[assay] == 'Hemi':
                        assay_dic[assay]['Seg'] += 1
                    elif row[assay] == 'Seg':
                        assay_dic[assay]['Seg'] += 1

                    elif row[assay] == 'Null':
                        assay_dic[assay]['Wildtype'] += 1
                    elif row[assay] == 'Wildtype':
                        assay_dic[assay]['Wildtype'] += 1

                    elif row[assay] == 'No Call':
                        assay_dic[assay]['No Call'] += 1

                    elif row[assay] == 'No Data':
                        assay_dic[assay]['Fail'] += 1

                    else:
                        continue

        for assay_id in self.headers:  # get the % data return key/value for the dictionary
            total = float(sum(assay_dic[assay_id].values()))
            total_calls = total - float(assay_dic[assay_id]['Fail'])
            percent_return = ((total_calls/total)*100)
            data_return = '%.2f' % percent_return

            assay_dic[assay_id]['% Data Return'] = float(data_return)

        return assay_dic  # return the dictionary
