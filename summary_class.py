class SummaryTable(object):
    """grab the reader data from kraken file and make a dictionary of the # of calls for each assay"""

    def __init__(self, kraken, headers):
        self.kraken = kraken
        self.headers = [i for i in headers if 'Zygosity' in i]

    def get_summary(self):
        """update counts of each assays calls and return the dictionary"""

        assay_dic = {}  # empty dictionary to append assays and calls to

        #  trait = G:G, seg = G:A, wildtype = A:A
        type_1 = ['DBSNP357223 Zygosity Call', 'DBSNP357253 Zygosity Call',
                  'DBSNP357256 Zygosity Call', '298535 Zygosity Call', '295699 Zygosity Call']

        #  trait = G:G, seg = G:T, wildtype = T:T
        type_2 = ['298512 Zygosity Call']

        #  trait = T:T, seg = G:T, wildtype = G:G
        type_3 = ['298516 Zygosity Call']

        #  homo = T:T, hemi = C:T, null = C:C
        type_4 = ['298518 Zygosity Call', 'DBSNP357255 Zygosity Call', '295698 Zygosity Call']

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

        for assay_id in self.headers:
            total = float(sum(assay_dic[assay_id].values()))
            total_calls = total - float(assay_dic[assay_id]['Fail'])
            percent_return = ((total_calls/total)*100)
            data_return = '%.2f' % percent_return

            assay_dic[assay_id]['% Data Return'] = float(data_return)

        return assay_dic
