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

        '''Types of assays and their respective allele calls as a Trait/Seg/Wildtype call'''

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

        """for each line in the file, count the trait/seg/null for the panel"""

        assay_num = int(len(self.dic))  # number of assays being run in the panel
        trait = 0
        seg = 0
        wildtype = 0
        no_call = 0
        no_data = 0

        for assay, allele in self.dic.items():

            if assay in type_1:  # Trait = G:G, Seg = G:A, Wildtype = A:A

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

            elif assay in type_2:  # Trait = G:G, Seg = G:T, Wildtype = T:T

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

            elif assay in type_3:  # Trait = T:T, Seg = G:T, Wildtype = G:G

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

            elif assay in type_4:  # Trait = T:T, Seg = C:T, Wildtype = C:C

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

            elif assay in type_5:  # Trait = T:T, Seg = T:A, Wildtype = A:A

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

            elif assay in type_6:  # Trait = A:A, Seg = G:A, Wildtype = G:G

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

            elif assay in type_7:  # Trait = C:C, Seg = C:A, Wildtype = A:A

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

            elif assay in type_8:  # Trait = G:G, Seg = G:C, Wildtype = C:C

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

            elif assay in type_9:  # Trait = C:C, Seg = C:T, Wildtype = T:T

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

            elif assay in type_10:  # Trait = A:A, Seg = T:A, Wildtype = 'T:T

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

            elif assay in type_11:  # Trait = A:A, Seg = C:A, Wildtype = C:C

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
            elif trait == (assay_num/2) and seg == (assay_num/2):
                tlm = 'SegT'
            elif wildtype == (assay_num/2) and seg == (assay_num/2):
                tlm = 'SegW'
            elif trait == (assay_num/2) and wildtype == (assay_num/2):
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
            elif trait == (assay_num-1) and wildtype == 1:
                tlm = 'BreakT'
            elif trait == 1 and wildtype == (assay_num-1):
                tlm = 'BreakW'
            elif trait == (assay_num/2) and wildtype == (assay_num/2):
                tlm = 'Undecided'
            else:
                tlm = 'Undecided'

            return tlm


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
        if self.call == 'Null':
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
