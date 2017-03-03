class TraitN9(object):
    """returns a call for the assay line in a kraken file if the N9 panel is present in the study"""

    def __init__(self, assay1, assay2, assay3, assay4):
        self.assay1 = assay1
        self.assay2 = assay2
        self.assay3 = assay3
        self.assay4 = assay4

    def n9_call(self):
        """determines the call on the N9 panel assays"""
        trait = 0
        seg = 0
        wildtype = 0
        no_call = 0
        no_data = 0

        #  get the number of homo/hemi/null calls for first assay in N9 panel -- 'DBSNP357223 Zygosity Call'
        if self.assay1 == 'G:G':
            trait += 1
        if self.assay1 == 'G:A':
            seg += 1
        if self.assay1 == 'A:A':
            wildtype += 1
        if self.assay1 == 'No Call':
            no_call += 1
        if self.assay1 == 'No Data':
            no_data += 1

        # get the number of homo/hemi/null calls for second assay in N9 panel -- 'DBSNP357253 Zygosity Call'
        if self.assay2 == 'G:G':
            trait += 1
        if self.assay2 == 'G:A':
            seg += 1
        if self.assay2 == 'A:A':
            wildtype += 1
        if self.assay2 == 'No Call':
            no_call += 1
        if self.assay2 == 'No Data':
            no_data += 1

            # get the number of homo/hemi/null calls for third assay in N9 panel -- 'DBSNP357255 Zygosity Call'
        if self.assay3 == 'T:T':
            trait += 1
        if self.assay3 == 'C:T':
            seg += 1
        if self.assay3 == 'C:C':
            wildtype += 1
        if self.assay3 == 'No Call':
            no_call += 1
        if self.assay3 == 'No Data':
            no_data += 1

        # get the number of homo/hemi/null calls for fourth assay in N9 panel -- 'DBSNP357256 Zygosity Call'
        if self.assay4 == 'G:G':
            trait += 1
        if self.assay4 == 'G:A':
            seg += 1
        if self.assay4 == 'A:A':
            wildtype += 1
        if self.assay4 == 'No Call':
            no_call += 1
        if self.assay4 == 'No Data':
            no_data += 1

        # get the call based on the number of homo/hemi/nulls per line
        if no_call > 0:
            nid = 'No Call'
        elif no_data > 0:
            nid = 'No Data'
        elif trait == 4:
            nid = 'Trait'
        elif seg == 4:
            nid = 'Seg'
        elif wildtype == 4:
            nid = 'Wildtype'
        elif trait in range(1, 4) and seg in range(1, 4) and wildtype == 0:
            nid = 'SegT'
        elif wildtype in range(1, 4) and seg in range(1, 4) and trait == 0:
            nid = 'SegW'
        elif trait == 2 and wildtype == 2:
            nid = 'Undecided'
        elif trait == 3 and wildtype == 1:
            nid = 'BreakT'
        elif trait == 1 and wildtype == 3:
            nid = 'BreakW'
        else:
            nid = 'Undecided'

        return nid


class TraitCastle(object):
    """returns a call for the assay line in a kraken file"""

    def __init__(self, assay1, assay2, assay3, assay4):
        self.assay1 = assay1
        self.assay2 = assay2
        self.assay3 = assay3
        self.assay4 = assay4

    def cas_call(self):
        """determines the call on the N9 panel assays"""
        trait = 0
        seg = 0
        wildtype = 0
        no_call = 0
        no_data = 0

        #  get the number of homo/hemi/null calls for first assay in N9 panel -- '298512 Zygosity Call'
        if self.assay1 == 'G:G':
            trait += 1
        if self.assay1 == 'G:T':
            seg += 1
        if self.assay1 == 'T:T':
            wildtype += 1
        if self.assay1 == 'No Call':
            no_call += 1
        if self.assay1 == 'No Data':
            no_data += 1

        # get the number of homo/hemi/null calls for second assay in N9 panel -- '298516 Zygosity Call'
        if self.assay2 == 'T:T':
            trait += 1
        if self.assay2 == 'G:T':
            seg += 1
        if self.assay2 == 'G:G':
            wildtype += 1
        if self.assay2 == 'No Call':
            no_call += 1
        if self.assay2 == 'No Data':
            no_data += 1

            # get the number of homo/hemi/null calls for third assay in N9 panel -- '298518 Zygosity Call'
        if self.assay3 == 'T:T':
            trait += 1
        if self.assay3 == 'C:T':
            seg += 1
        if self.assay3 == 'C:C':
            wildtype += 1
        if self.assay3 == 'No Call':
            no_call += 1
        if self.assay3 == 'No Data':
            no_data += 1

        # get the number of homo/hemi/null calls for fourth assay in N9 panel -- '298535 Zygosity Call'
        if self.assay4 == 'G:G':
            trait += 1
        if self.assay4 == 'G:A':
            seg += 1
        if self.assay4 == 'A:A':
            wildtype += 1
        if self.assay4 == 'No Call':
            no_call += 1
        if self.assay4 == 'No Data':
            no_data += 1

        # get the call based on the number of homo/hemi/nulls per line
        if no_call > 0:
            casid = 'No Call'
        elif no_data > 0:
            casid = 'No Data'
        elif trait == 4:
            casid = 'Trait'
        elif seg == 4:
            casid = 'Seg'
        elif wildtype == 4:
            casid = 'Wildtype'
        elif trait in range(1, 4) and seg in range(1, 4) and wildtype == 0:
            casid = 'SegT'
        elif wildtype in range(1, 4) and seg in range(1, 4) and trait == 0:
            casid = 'SegW'
        elif trait == 2 and wildtype == 2:
            casid = 'Undecided'
        elif trait == 3 and wildtype == 1:
            casid = 'BreakT'
        elif trait == 1 and wildtype == 3:
            casid = 'BreakW'
        else:
            casid = 'Undecided'

        return casid


class TraitCR(object):
    """returns a call for the assay line in a kraken file"""

    def __init__(self, assay1, assay2):
        self.assay1 = assay1
        self.assay2 = assay2

    def cr_call(self):
        """determines the call on the CR panel assays"""
        homo = 0
        hemi = 0
        null = 0
        no_call = 0
        no_data = 0

        #  get the number of homo/hemi/null calls for first assay in N9 panel -- 'DBSNP357223 Zygosity Call'
        if self.assay1 == 'Homo':
            homo += 1
        if self.assay1 == 'Hemi':
            hemi += 1
        if self.assay1 == 'Null':
            null += 1
        if self.assay1 == 'No Call':
            no_call += 1
        if self.assay1 == 'No Data':
            no_data += 1

        # get the number of homo/hemi/null calls for second assay in N9 panel -- 'DBSNP357253 Zygosity Call'
        if self.assay2 == 'Homo':
            homo += 1
        if self.assay2 == 'Hemi':
            hemi += 1
        if self.assay2 == 'Null':
            null += 1
        if self.assay2 == 'No Call':
            no_call += 1
        if self.assay2 == 'No Data':
            no_data += 1

        # get the call based on the number of homo/hemi/nulls per line
        if no_call > 0:
            crid = 'No Call'
        elif no_data > 0:
            crid = 'No Data'
        elif homo == 2:
            crid = 'Trait'
        elif hemi == 2:
            crid = 'Seg'
        elif null == 2:
            crid = 'Wildtype'
        elif homo == 1 and hemi == 1:
            crid = 'SegT'
        elif null == 1 and hemi == 1:
            crid = 'SegW'
        elif homo == 1 and null == 1:
            crid = 'Undecided'
        else:
            crid = 'Undecided'

        return crid


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
