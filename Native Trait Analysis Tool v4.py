# -*- coding: utf-8 -*-
"""
This is version 4 of the TLM tool - author Brian James (u590135)

This version is written in python 3.5
"""

import csv
import pandas as pd
from pandas import ExcelWriter
import glob
import os
import time

"""dictionary of all assays in production and their respective call of Trait/Seg/Wildtype for each allele.
Any SNP addition to the panel needs to have its name added to the type list corresponding to the allele
call for that assay"""
assay_types = {
    # type1 : Trait = G:G, Seg = G:A, Wildtype = A:A
    'type1': [
        'DBSNP357223 Zygosity Call',
        'DBSNP357253 Zygosity Call',
        'DBSNP357256 Zygosity Call',
        '298535 Zygosity Call',
        '295699 Zygosity Call',
        '303605 Zygosity Call',
        'DBSNP357252 Zygosity Call',
        'N13-7444437 Zygosity Call'
    ],

    # type2 : Trait = G:G, Seg = G:T, Wildtype = T:T
    'type2': [
        '298512 Zygosity Call'
    ],

    # type3 : Trait = T:T, Seg = G:T, Wildtype = G:G
    'type3': [
        '298516 Zygosity Call',
        '303806 Zygosity Call'
    ],

    # type4 : Trait = T:T, Seg = C:T, Wildtype = C:C
    'type4': [
        '298518 Zygosity Call',
        'DBSNP357255 Zygosity Call',
        '295698 Zygosity Call',
        '302267 Zygosity Call',
        '281051 Zygosity Call',
        '220074 Zygosity Call',
        '320752 Zygosity Call'
    ],

    # type5 : Trait = T:T, Seg = T:A, Wildtype = A:A
    'type5': [
        '15473 Zygosity Call'
    ],

    # type6 : Trait = A:A, Seg = G:A, Wildtype = G:G
    'type6': [
        '62768 Zygosity Call',
        '295723 Zygosity Call',
        'N13-7444161 Zygosity Call',
        '335227 Zygosity Call'
    ],

    # type7 : Trait = C:C, Seg = C:A, Wildtype = A:A
    'type7': [
        '331689 Zygosity Call',
        '303572 Zygosity Call',
        'N13-7444297 Zygosity Call'
    ],

    # type8 : Trait = G:G, Seg = G:C, Wildtype = C:C
    'type8': [
        'DBSNP357226 Zygosity Call',
        '15925 Zygosity Call'
    ],

    # type9 : Trait = C:C, Seg = C:T, Wildtype = T:T
    'type9': [
        'DBSNP357233 Zygosity Call',
        'DBSNP357257 Zygosity Call',
        '332445 Zygosity Call',
        '320634 Zygosity Call'
    ],

    # type10 : Trait = A:A, Seg = T:A, Wildtype = 'T:T
    'type10': [
        '320826 Zygosity Call',
        '106265 Zygosity Call'
    ],

    # type11 : Trait = A:A, Seg = C:A, Wildtype = C:C
    'type11': [
        '328092 Zygosity Call'
    ]
}


def indexAssays(panel, headers):
    """read through Kraken file and return an index of where that panel is located in file.
    If there are additions to the panel, they must be added here to their respective variable"""

    panel_assays = {
        'N9_Panel': ['DBSNP357223 Zygosity Call',
                     'DBSNP357226 Zygosity Call',
                     'DBSNP357233 Zygosity Call',
                     'DBSNP357252 Zygosity Call',
                     'DBSNP357253 Zygosity Call',
                     'DBSNP357255 Zygosity Call',
                     'DBSNP357256 Zygosity Call',
                     'DBSNP357257 Zygosity Call'],

        'CastlePanel': ['298512 Zygosity Call',
                        '298516 Zygosity Call',
                        '298518 Zygosity Call',
                        '298535 Zygosity Call'],

        'CR_Panel': ['CRM1 Zygosity Call',
                     'CRM2 Zygosity Call'],

        'RLM7_Spring': ['295698 Zygosity Call',
                        '295699 Zygosity Call'],

        'RLM7_Winter': ['15473 Zygosity Call',
                        '62768 Zygosity Call',
                        '295723 Zygosity Call'],

        'BL10_Panel': ['302267 Zygosity Call',
                       '331689 Zygosity Call',
                       '281051 Zygosity Call',
                       '303572 Zygosity Call',
                       '303605 Zygosity Call',
                       '303806 Zygosity Call'],

        'N13_Panel': ['15925 Zygosity Call',
                      '106265 Zygosity Call',
                      '220074 Zygosity Call',
                      '320634 Zygosity Call',
                      '320752 Zygosity Call',
                      '320826 Zygosity Call',
                      '328092 Zygosity Call',
                      '332445 Zygosity Call',
                      '335227 Zygosity Call',
                      'N13-7444161 Zygosity Call',
                      'N13-7444297 Zygosity Call',
                      'N13-7444437 Zygosity Call'],

        'FAE_Panel': ['FAE 1-1 A Zygosity Call',
                      'FAE 1-2 C Zygosity Call']
    }

    panel_index = []

    for i in panel_assays[panel]:
        # noinspection PyBroadException
        try:
            panel_index.append(headers.index(i))  # get index of N9 panel (if exists)
        except:
            continue
    return panel_index


class SummaryTable(object):
    """grab the reader data from kraken file and make a dictionary of the # of calls for each assay"""

    def __init__(self, kraken, headers):
        self.kraken = kraken
        # only grab headers that have 'zygosity' in name
        self.headers = [i for i in headers if 'Zygosity' in i]

    def get_summary(self):
        """update counts of each assays calls and return the dictionary"""
        global assay_types

        assay_dic = {}  # empty dictionary to append assays and calls to

        #  append all assays in kraken file to the dictionary with 0 as blank call to begin with
        for i in self.headers:
            assay_dic[i] = {'Trait': 0, 'Seg': 0, 'Wildtype': 0, 'No Call': 0, 'Fail': 0}

        # read the kraken study line by line and start counting the homo/hemi/null values in the dictionary for
        #  each assay
        for row in self.kraken:

            for assay in self.headers:

                if assay in assay_types['type1']:
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

                elif assay in assay_types['type2']:
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

                elif assay in assay_types['type3']:
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

                elif assay in assay_types['type4']:
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

                elif assay in assay_types['type5']:
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

                elif assay in assay_types['type6']:
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

                elif assay in assay_types['type7']:
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

                elif assay in assay_types['type8']:
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

                elif assay in assay_types['type9']:
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

                elif assay in assay_types['type10']:
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

                elif assay in assay_types['type11']:
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
            percent_return = ((total_calls / total) * 100)
            data_return = '%.2f' % percent_return

            assay_dic[assay_id]['% Data Return'] = float(data_return)

        return assay_dic  # return the dictionary


class TraitLinkedMarker(object):
    """create a call for any of the Trait Linked Marker panels for the MAS group"""

    def __init__(self, assay_dict):
        self.dic = assay_dict
        assay_count = 1
        for k, v in assay_dict.items():
            assay = 'assay' + str(assay_count)
            setattr(self, assay, v)
            assay_count += 1

    def get_tlm_call(self):
        """count the trait/seg/wildtype ratios for any panel initiated in class and return a corresponding call"""
        global assay_types

        """for each line in the file, count the trait/seg/null for the panel"""

        assay_num = int(len(self.dic))  # number of assays being run in the panel
        trait = 0
        seg = 0
        wildtype = 0
        no_call = 0
        no_data = 0

        for assay, allele in self.dic.items():

            if assay in assay_types['type1']:

                if allele == 'G:G':
                    trait += 1
                elif allele == 'G:A':
                    seg += 1
                elif allele == 'A:G':
                    seg += 1
                elif allele == 'A:A':
                    wildtype += 1
                elif allele == 'No Call':
                    no_call += 1
                elif allele == 'No Data':
                    no_data += 1

            elif assay in assay_types['type2']:

                if allele == 'G:G':
                    trait += 1
                elif allele == 'G:T':
                    seg += 1
                elif allele == 'T:G':
                    seg += 1
                elif allele == 'T:T':
                    wildtype += 1
                elif allele == 'No Call':
                    no_call += 1
                elif allele == 'No Data':
                    no_data += 1

            elif assay in assay_types['type3']:

                if allele == 'T:T':
                    trait += 1
                elif allele == 'G:T':
                    seg += 1
                elif allele == 'T:G':
                    seg += 1
                elif allele == 'G:G':
                    wildtype += 1
                elif allele == 'No Call':
                    no_call += 1
                elif allele == 'No Data':
                    no_data += 1

            elif assay in assay_types['type4']:

                if allele == 'T:T':
                    trait += 1
                elif allele == 'C:T':
                    seg += 1
                elif allele == 'C:T':
                    seg += 1
                elif allele == 'C:C':
                    wildtype += 1
                elif allele == 'No Call':
                    no_call += 1
                elif allele == 'No Data':
                    no_data += 1

            elif assay in assay_types['type5']:

                if allele == 'T:T':
                    trait += 1
                elif allele == 'A:T':
                    seg += 1
                elif allele == 'T:A':
                    seg += 1
                elif allele == 'A:A':
                    wildtype += 1
                elif allele == 'No Call':
                    no_call += 1
                elif allele == 'No Data':
                    no_data += 1

            elif assay in assay_types['type6']:

                if allele == 'A:A':
                    trait += 1
                elif allele == 'G:A':
                    seg += 1
                elif allele == 'A:G':
                    seg += 1
                elif allele == 'G:G':
                    wildtype += 1
                elif allele == 'No Call':
                    no_call += 1
                elif allele == 'No Data':
                    no_data += 1

            elif assay in assay_types['type7']:

                if allele == 'C:C':
                    trait += 1
                elif allele == 'C:A':
                    seg += 1
                elif allele == 'A:C':
                    seg += 1
                elif allele == 'A:A':
                    wildtype += 1
                elif allele == 'No Call':
                    no_call += 1
                elif allele == 'No Data':
                    no_data += 1

            elif assay in assay_types['type8']:

                if allele == 'G:G':
                    trait += 1
                elif allele == 'G:C':
                    seg += 1
                elif allele == 'C:G':
                    seg += 1
                elif allele == 'C:C':
                    wildtype += 1
                elif allele == 'No Call':
                    no_call += 1
                elif allele == 'No Data':
                    no_data += 1

            elif assay in assay_types['type9']:

                if allele == 'C:C':
                    trait += 1
                elif allele == 'C:T':
                    seg += 1
                elif allele == 'T:C':
                    seg += 1
                elif allele == 'T:T':
                    wildtype += 1
                elif allele == 'No Call':
                    no_call += 1
                elif allele == 'No Data':
                    no_data += 1

            elif assay in assay_types['type10']:

                if allele == 'A:A':
                    trait += 1
                elif allele == 'T:A':
                    seg += 1
                elif allele == 'A:T':
                    seg += 1
                elif allele == 'T:T':
                    wildtype += 1
                elif allele == 'No Call':
                    no_call += 1
                elif allele == 'No Data':
                    no_data += 1

            elif assay in assay_types['type11']:

                if allele == 'A:A':
                    trait += 1
                elif allele == 'C:A':
                    seg += 1
                elif allele == 'A:C':
                    seg += 1
                elif allele == 'C:C':
                    wildtype += 1
                elif allele == 'No Call':
                    no_call += 1
                elif allele == 'No Data':
                    no_data += 1

            else:

                if allele == 'Trait':
                    trait += 1
                elif allele == 'Homo':
                    trait += 1
                elif allele == 'Seg':
                    seg += 1
                elif allele == 'Hemi':
                    seg += 1
                elif allele == 'Wildtype':
                    wildtype += 1
                elif allele == 'Null':
                    wildtype += 1
                elif allele == 'No Call':
                    no_call += 1
                elif allele == 'No Data':
                    no_data += 1

        """Using the ratios from the counting above, come up with a returnable call as Trait/Seg/Wildtype"""

        if assay_num == 2:  # if only 2 assays in TLM panel

            if no_data > 0:
                tlm = 'No Data'
            elif no_call > 0:
                tlm = 'No Call'
            elif trait == assay_num:
                tlm = 'Trait'
            elif seg == assay_num:
                tlm = 'Seg'
            elif wildtype == assay_num:
                tlm = 'Wildtype'
            elif trait == (assay_num / 2) and seg == (assay_num / 2):
                tlm = 'SegT'
            elif wildtype == (assay_num / 2) and seg == (assay_num / 2):
                tlm = 'SegW'
            elif trait == (assay_num / 2) and wildtype == (assay_num / 2):
                tlm = 'Undecided'
            else:
                tlm = 'Undecided'

            return tlm

        else:  # if the number of assays in the TLM panel is greater than 2

            if no_data > 0:
                tlm = 'No Data'
            elif no_call > 0:
                tlm = 'No Call'
            elif trait == assay_num:
                tlm = 'Trait'
            elif seg == assay_num:
                tlm = 'Seg'
            elif wildtype == assay_num:
                tlm = 'Wildtype'
            elif trait in range(1, assay_num) and seg in range(1, assay_num) and wildtype == 0:
                tlm = 'SegT'
            elif wildtype in range(1, assay_num) and seg in range(1, assay_num) and trait == 0:
                tlm = 'SegW'
            elif trait == (assay_num - 1) and wildtype == 1:
                tlm = 'BreakT'
            elif trait == 1 and wildtype == (assay_num - 1):
                tlm = 'BreakW'
            elif trait == (assay_num / 2) and wildtype == (assay_num / 2):
                tlm = 'Undecided'
            else:
                tlm = 'Undecided'

            return tlm

    def get_fae(self):
        """Read the FAE 1-1A and FAE 1-2C markers and return a summarized call"""

        trait = 0
        seg = 0
        wildtype = 0
        no_call = 0
        no_data = 0

        for assay, allele in self.dic.items():

            if assay == 'FAE 1-1 A Zygosity Call':
                if allele == 'Trait':
                    trait += 1
                elif allele == 'Homo':
                    trait += 1
                elif allele == 'Seg':
                    seg += 1
                elif allele == 'Hemi':
                    seg += 1
                elif allele == 'Wildtype':
                    wildtype += 1
                elif allele == 'Null':
                    wildtype += 1
                elif allele == 'No Call':
                    no_call += 1
                elif allele == 'No Data':
                    no_data += 1

            elif assay == 'FAE 1-2 C Zygosity Call':
                if allele == 'Trait':
                    trait += 1
                elif allele == 'Homo':
                    trait += 1
                elif allele == 'Seg':
                    seg += 1
                elif allele == 'Hemi':
                    seg += 1
                elif allele == 'Wildtype':
                    wildtype += 1
                elif allele == 'Null':
                    wildtype += 1
                elif allele == 'No Call':
                    no_call += 1
                elif allele == 'No Data':
                    no_data += 1

            '''return the call based on the following criteria'''

            if no_data > 0:
                return 'No Data'
            elif no_call > 0:
                return 'No Call'
            elif wildtype > 0:
                return 'Wildtype'
            elif seg in range(1, 3) and wildtype == 0:
                return 'Seg'
            elif trait == 2:
                return 'Trait'


class CallConversion(object):
    """converts a homo/hemi/null call into the trait/hemi/wildtype or plus/minus calls for Variety book"""

    def __init__(self, call):
        self.call = call

    def zygo_conv(self):
        if self.call == 'Homo':
            return 'Trait'
        if self.call == 'Hemi':
            return 'Seg'
        if self.call == 'Null':
            return 'Wildtype'
        else:
            return self.call

    def plus_minus(self):
        if self.call == 'Homo':
            return 'Plus'
        if self.call == 'Trait':
            return 'Plus'
        if self.call == 'Null':
            return 'Minus'
        if self.call == 'Wildtype':
            return 'Minus'
        else:
            return self.call

    def cyto_call(self):
        if self.call == 'Trait':
            return 'A-Cyto'
        if self.call == 'Wildtype':
            return 'B-Cyto'
        else:
            return self.call


# ------------------------------------- finding the files in the folder ----------------------------
'''reading the kraken csv file with pandas and converting the data set into individual dataframes'''
#fhand = 'testfile.csv'
path = r'C:\Users\U590135\.spyder-py3\projects\TLM'
xlextension = 'xlsx'
extension = 'csv'
os.chdir(path)
xlsx_result = [xl for xl in glob.glob('*.{}'.format(xlextension))]

# convert the xlsx files into csv files
xl_files = []
for xlsx_files in xlsx_result:
    if '_report' in xlsx_files:
        continue
    else:
        xl_files.append(xlsx_files)

for i in xl_files:
    df_temp = pd.read_excel(i, index=False, encoding='utf-8')
    get_name = i.find('.')
    csv_name = i[:get_name] + '.csv'
    df_temp.to_csv(csv_name, index=False, encoding='utf-8')
    #os.remove(i)

#  get a list of all the csv files to run the script on
csv_result = [i for i in glob.glob('*.{}'.format(extension))]

completed_path = r'C:\Users\U590135\.spyder-py3\projects\TLM\completed'

results_files = []
for i in csv_result:
    if 'temp_calls' in i:
        continue
    else:
        results_files.append(i)

# ------------------------------------- Summary Table Creation -------------------------------------

for fhand in results_files:
    print('Analyzing ' + fhand + '...')

    with open(fhand, newline='') as kraken:
        reader = csv.DictReader(kraken)
        kraken_headers = reader.fieldnames
                
        summary_assays = [i for i in kraken_headers if 'Zygosity' in i]  # select only assays for summary
        summary_init = SummaryTable(reader, kraken_headers)
        summary_dict = summary_init.get_summary()
        sum_df_initial = pd.DataFrame.from_dict(summary_dict)
        sum_df_assays = sum_df_initial[summary_assays]
        summary_df = sum_df_assays.reindex(['Trait', 'Seg', 'Wildtype', 'No Call',
                                            'Fail', '% Data Return'])
        
    # --------------------------------------- Isolating the Panels -------------------------------------
    
        # Notes - the next few sections will read the kraken file and isolate each panel (if present) by
        # trying to index the assays.  If it is in the file, the index will tell us which columns belong to
        # which panel so we can isolate the information to pass into TLM class and make a summary
        # column in reader.
    
    with open(fhand, newline='') as kraken:
        reader = csv.DictReader(kraken)
        kraken_headers = reader.fieldnames
    
        panel_summaries = ['ACM N9 Summary', 'Castle Summary', 'CR Mendel Summary', 'RLM 7 Spring Summary',
                           'RLM 7 Winter Summary', 'BL10 Summary', 'ACM N13 Summary', 'FAE Summary Call']
    
        panel_list = []  # empty list where we can store which panels were run - cxzfor final merge sequence
        header_list = kraken_headers

        for summs in panel_summaries:
            header_list.append(summs)  # add the possible panel summaries to the headers
    
        with open('temp_calls.csv', 'w', newline='') as resultfile:
            writer = csv.DictWriter(resultfile, fieldnames=header_list)
            writer.writeheader()
    
            # try and index all of the panels to see if they are in the kraken analysis
    
            n9_index = indexAssays('N9_Panel', kraken_headers)
            n9_len = int(len(n9_index))
            if n9_len > 1:
                panel_list.append('N9_Panel')
    
            cas_index = indexAssays('CastlePanel', kraken_headers)
            cas_len = int(len(cas_index))
            if cas_len > 1:
                panel_list.append('CastlePanel')
    
            cr_index = indexAssays('CR_Panel', kraken_headers)
            cr_len = int(len(cr_index))
            if cr_len > 1:
                panel_list.append('CR_Panel')
    
            rlm7s_index = indexAssays('RLM7_Spring', kraken_headers)
            rlm7s_len = int(len(rlm7s_index))
            if rlm7s_len > 1:
                panel_list.append('RLM7_Spring')
    
            rlm7w_index = indexAssays('RLM7_Winter', kraken_headers)
            rlm7w_len = int(len(rlm7w_index))
            if rlm7w_len > 1:
                panel_list.append('RLM7_Winter')
    
            bl10_index = indexAssays('BL10_Panel', kraken_headers)
            bl10_len = int(len(bl10_index))
            if bl10_len > 1:
                panel_list.append('BL10_Panel')
    
            n13_index = indexAssays('N13_Panel', kraken_headers)
            n13_len = int(len(n13_index))
            if n13_len > 1:
                panel_list.append('N13_Panel')
    
            fae_index = indexAssays('FAE_Panel', kraken_headers)
            fae_len = int(len(fae_index))
            if fae_len > 1:
                panel_list.append('FAE_Panel')
            
            for row in reader:
                remove_controls = ['H7', 'H07', 'H8', 'H08', 'H9', 'H09', 'H10', 'H11', 'H12']
                if row['Well'] in remove_controls:  # remove controls from report
                    continue
                else:    
                    row['Pedigree'] = row['Pedigree'].encode('utf-8')  # encode the pedigree column as utf-8
    
                # ------------------------  N9 Panel ----------------------------------
                if n9_len > 1:
                    ''' N9 is present so make a dictionary out of the panel, pass into TLM to make a summary column'''
                    n9_index.sort()  # make sure column names are in order
                    dict_n9 = {}
    
                    for i in n9_index:
                        dict_n9[kraken_headers[i]] = row[kraken_headers[i]]
    
                    n9 = TraitLinkedMarker(dict_n9)
                    call_n9 = n9.get_tlm_call()
                    row['ACM N9 Summary'] = call_n9
    
                # --------------------------- Castle Panel --------------------------------
                if cas_len > 1:
                    ''' Castle is present so make a dictionary out of the panel, pass into TLM to make a summary column'''
                    cas_index.sort()  # make sure column names are in order
                    dict_cas = {}
    
                    for i in cas_index:
                        dict_cas[kraken_headers[i]] = row[kraken_headers[i]]
    
                    castle = TraitLinkedMarker(dict_cas)
                    call_cas = castle.get_tlm_call()
                    row['Castle Summary'] = call_cas
    
                # ------------------------------ CR Panel ---------------------------------
                if cr_len > 1:
                    ''' CR is present so make a dictionary out of the panel, pass into TLM to make a summary column'''
                    cr_index.sort()  # make sure column names are in order
                    dict_cr = {}
    
                    for i in cr_index:
                        dict_cr[kraken_headers[i]] = row[kraken_headers[i]]
    
                    cr = TraitLinkedMarker(dict_cr)
                    call_cr = cr.get_tlm_call()
                    row['CR Mendel Summary'] = call_cr
    
                # ---------------------------- RLM7 Spring Panel --------------------------
                if rlm7s_len > 1:
                    ''' RLM7 Spring is present so convert panel into dictionary, pass into TLM to make a summary column'''
                    rlm7s_index.sort()  # make sure column names are in order
                    dict_rlm7s = {}
    
                    for i in rlm7s_index:
                        dict_rlm7s[kraken_headers[i]] = row[kraken_headers[i]]
    
                    rlm7s = TraitLinkedMarker(dict_rlm7s)
                    call_rlm7s = rlm7s.get_tlm_call()
                    row['RLM 7 Spring Summary'] = call_rlm7s
    
                # ------------------------------ RLM7 Winter Panel ------------------------
                if rlm7w_len > 1:
                    ''' RLM7 Winter is present so convert panel into dictionary, pass into TLM to make a summary column'''
                    rlm7w_index.sort()  # make sure column names are in order
                    dict_rlm7w = {}
    
                    for i in rlm7w_index:
                        dict_rlm7w[kraken_headers[i]] = row[kraken_headers[i]]
    
                    rlm7w = TraitLinkedMarker(dict_rlm7w)
                    call_rlm7w = rlm7w.get_tlm_call()
                    row['RLM 7 Winter Summary'] = call_rlm7w
    
                # ------------------------------- BL10 Panel ------------------------------
                if bl10_len > 1:
                    ''' BL10 is present so convert panel into dictionary, pass into TLM to make a summary column'''
                    bl10_index.sort()  # make sure column names are in order
                    dict_bl10 = {}
    
                    for i in bl10_index:
                        dict_bl10[kraken_headers[i]] = row[kraken_headers[i]]
    
                    bl10 = TraitLinkedMarker(dict_bl10)
                    call_bl10 = bl10.get_tlm_call()
                    row['BL10 Summary'] = call_bl10
    
                # ------------------------------- N13 Panel -------------------------------
                if n13_len > 1:
                    ''' N13 is present so convert panel into dictionary, pass into TLM to make a summary column'''
                    n13_index.sort()  # make sure column names are in order
                    dict_n13 = {}
    
                    for i in n13_index:
                        dict_n13[kraken_headers[i]] = row[kraken_headers[i]]
    
                    n13 = TraitLinkedMarker(dict_n13)
                    call_n13 = n13.get_tlm_call()
                    row['ACM N13 Summary'] = call_n13
    
                # ------------------------------- FAE Panel -------------------------------
                if fae_len > 0:
                    ''' FAE is present so convert panel into dictionary, pass into TLM to make a summary column'''
                    fae_index.sort()  # make sure column names are in order
                    dict_fae = {}
    
                    for i in fae_index:
                        dict_fae[kraken_headers[i]] = row[kraken_headers[i]]
    
                    fae = TraitLinkedMarker(dict_fae)
                    call_fae = fae.get_fae()
                    row['FAE Summary Call'] = call_fae
    
                writer.writerow(row)
    
    # -----------format the file so summary columns are after assay panels, and replace cyto/GT200/lepr3 -------------------
    
    tlm_df = pd.read_csv('temp_calls.csv', encoding='UTF-8')
    # get rid of the b' prefix in pedigree column
    tlm_df['Pedigree'] = tlm_df['Pedigree'].map(lambda x: x.lstrip("b'").rstrip("'"))

    
    #  replace the column names and contents for the following assays if present in the data file
    if 'Cyto Zygosity Call' in tlm_df.columns:
        tlm_df = tlm_df.replace(to_replace={'Cyto Zygosity Call': {'Trait': 'A-Cyto'}})
        tlm_df = tlm_df.replace(to_replace={'Cyto Zygosity Call': {'Wildtype': 'B-Cyto'}})
        tlm_df.rename(columns={'Cyto Zygosity Call': 'Cyto Call'}, inplace=True)
    
    if 'GT 200 Zygosity Call' in tlm_df.columns:
        tlm_df = tlm_df.replace(to_replace={'GT 200 Zygosity Call': {'Homo': 'Plus'}})
        tlm_df = tlm_df.replace(to_replace={'GT 200 Zygosity Call': {'Null': 'Minus'}})
        tlm_df.rename(columns={'GT 200 Zygosity Call': 'GT200 Call'}, inplace=True)
    
    if 'GT200 Zygosity Call' in tlm_df.columns:
        tlm_df = tlm_df.replace(to_replace={'GT200 Zygosity Call': {'Homo': 'Plus'}})
        tlm_df = tlm_df.replace(to_replace={'GT200 Zygosity Call': {'Null': 'Minus'}})
        tlm_df.rename(columns={'GT200 Zygosity Call': 'GT200 Call'}, inplace=True)
    
    if 'LepR3-R3S3A Zygosity Call' in tlm_df.columns:
        tlm_df = tlm_df.replace(to_replace={'LepR3-R3S3A Zygosity Call': {'Homo': 'Trait'}})
        tlm_df = tlm_df.replace(to_replace={'LepR3-R3S3A Zygosity Call': {'Hemi': 'Seg'}})
        tlm_df = tlm_df.replace(to_replace={'LepR3-R3S3A Zygosity Call': {'Null': 'Wildtype'}})
    
    if 'BGSP-02 Zygosity Call' in tlm_df.columns:
        tlm_df = tlm_df.replace(to_replace={'BGSP-02 Zygosity Call': {'Trait': 'Plus'}})
        tlm_df = tlm_df.replace(to_replace={'BGSP-02 Zygosity Call': {'Wildtype': 'Minus'}})
        tlm_df.rename(columns={'BGSP-02 Zygosity Call': 'BGSP-02'}, inplace=True)
    
    # reformat the file and write the TLM summary report and summary table to excel
    
    krak_info = ['Box', 'Well', 'Project', 'Pedigree', 'Source ID', 'Geno_Id', 'RowId', 'Loc Seq#']
    
    df_initial = tlm_df[krak_info]  # sample and well info dataframe only (no molecular data)
    
    # merge initial with N9 dataframe if panel is in the kraken file
    
    if 'N9_Panel' in panel_list:
        n9_headers = krak_info
        for i in n9_index:
            n9_headers.append(kraken_headers[i])
        n9_headers.append('ACM N9 Summary')
        n9_df = tlm_df[n9_headers]
        merge1 = pd.merge(left=df_initial, right=n9_df, how='left')
    else:
        tlm_df.drop(['ACM N9 Summary'], axis=1, inplace=True, errors='ignore')
        merge1 = df_initial
    
    # merge that with N13 dataframe if panel is in the kraken file
    
    if 'N13_Panel' in panel_list:
        n13_headers = krak_info
        for i in n13_index:
            n13_headers.append(kraken_headers[i])
        n13_headers.append('ACM N13 Summary')
        n13_df = tlm_df[n13_headers]
        merge2 = pd.merge(left=merge1, right=n13_df, how='left')
    else:
        tlm_df.drop(['ACM N13 Summary'], axis=1, inplace=True, errors='ignore')
        merge2 = merge1
    
    # merge that with castle dataframe if panel is in the kraken file
    
    if 'CastlePanel' in panel_list:
        cas_headers = krak_info
        for i in cas_index:
            cas_headers.append(kraken_headers[i])
        cas_headers.append('Castle Summary')
        cas_df = tlm_df[cas_headers]
        merge3 = pd.merge(left=merge2, right=cas_df, how='left')
    else:
        tlm_df.drop(['Castle Summary'], axis=1, inplace=True, errors='ignore')
        merge3 = merge2
    
    # merge that with clubroot dataframe if panel is in the kraken file
    
    if 'CR_Panel' in panel_list:
        cr_headers = krak_info
        for i in cr_index:
            cr_headers.append(kraken_headers[i])
        cr_headers.append('CR Mendel Summary')
        cr_df = tlm_df[cr_headers]
        merge4 = pd.merge(left=merge3, right=cr_df, how='left')
    else:
        tlm_df.drop(['CR Mendel Summary'], axis=1, inplace=True, errors='ignore')
        merge4 = merge3
    
    # merge that with RLM7 spring dataframe if panel is in the kraken file
    
    if 'RLM7_Spring' in panel_list:
        rlm7s_headers = krak_info
        for i in rlm7s_index:
            rlm7s_headers.append(kraken_headers[i])
        rlm7s_headers.append('RLM 7 Spring Summary')
        rlm7s_df = tlm_df[rlm7s_headers]
        merge5 = pd.merge(left=merge4, right=rlm7s_df, how='left')
    else:
        tlm_df.drop(['RLM 7 Spring Summary'], axis=1, inplace=True, errors='ignore')
        merge5 = merge4
    
    # merge that with RLM7 winter dataframe if that panel is in the kraken file
    
    if 'RLM7_Winter' in panel_list:
        rlm7w_headers = krak_info
        for i in rlm7w_index:
            rlm7w_headers.append(kraken_headers[i])
        rlm7w_headers.append('RLM 7 Winter Summary')
        rlm7w_df = tlm_df[rlm7w_headers]
        merge6 = pd.merge(left=merge5, right=rlm7w_df, how='left')
    else:
        tlm_df.drop(['RLM 7 Winter Summary'], axis=1, inplace=True, errors='ignore')
        merge6 = merge5
    
    # merge that with BL10 dataframe if that panel is in the kraken file
    
    if 'BL10_Panel' in panel_list:
        bl10_headers = krak_info
        for i in bl10_index:
            bl10_headers.append(kraken_headers[i])
        bl10_headers.append('BL10 Summary')
        bl10_df = tlm_df[bl10_headers]
        merge7 = pd.merge(left=merge6, right=bl10_df, how='left')
    else:
        tlm_df.drop(['BL10 Summary'], axis=1, inplace=True, errors='ignore')
        merge7 = merge6
    
    # merge that with FAE dataframe if that panel is in the kraken file
    
    if 'FAE_Panel' in panel_list:
        fae_headers = krak_info
        for i in fae_index:
            fae_headers.append(kraken_headers[i])
        fae_headers.append('FAE Summary Call')
        fae_df = tlm_df[fae_headers]
        merge8 = pd.merge(left=merge7, right=fae_df, how='left')
    else:
        tlm_df.drop(['FAE Summary Call'], axis=1, inplace=True, errors='ignore')
        merge8 = merge7
    
    # merge with final dataframe and write the file (gets any assays not in panel tacked on to end of kraken report)
    
    merge_final = pd.merge(left=merge8, right=tlm_df, how='left')
    
    find_the_dot = fhand.find('.')
    final_name = fhand[:find_the_dot] + '_Report.xlsx'
    results_writer = ExcelWriter(final_name)
    
    os.chdir(completed_path)
    merge_final.to_excel(results_writer, 'Report', index=False)
    summary_df.to_excel(results_writer, 'Summary Table')
    results_writer.save()
    
    os.chdir(path) 
    os.remove('temp_calls.csv')
    os.remove(fhand)
    
    print('Successfully completed {} ...\n'.format(fhand))
    
    time.sleep(2)
                
                
            
            
            
            
     
            
            
            