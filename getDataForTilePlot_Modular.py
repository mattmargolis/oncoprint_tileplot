#!/usr/bin/env python2
import sys
import os, logging, csv, collections
import argparse
from researchmaster.report import *

'''
Type is a digit code depending on the variant
    SV:
        4 = (steelblue4) nonsense or frameshift or splice
        1 = (orangered) all other SV
    CN:
        2 = (deep pink) amplification
        3 = (skyblue) deletion
    RE:
        4 = (steelblue4) truncations
        5 = (khaki) all other RE

    Multiple:
        6 = (red) amp + sv (2+1)
        7 = (blue) del + truncs (3+4)
        8 = (yellow) amp + RE (2+5)
        9 = (grey70) other multiple
'''

class OncoPrinter():
    def __init__(self, varDict, specimens, queryGenes,
                 GroupDict={}, GroupGeneList=list(), GroupsToSkip=list(), TopTen=False, MSI_only=False,
                PathologyList=list(), HRD_only=False, exclude_VUS=True, special_type_list=list()):
        self.varDict=varDict
        self.specimens=specimens

        self.queryGenes = queryGenes
        #queryGenes are list of Genes to plot for tileplotter
        self.noMutSamples=set()
        self.mutSamples=set()
        self.GroupDict=GroupDict
        self.GroupGeneList=GroupGeneList
        #GroupGeneList is the list of genes to skip because are included in groups
        self.GroupsToSkip=GroupsToSkip
        self.TopTen=TopTen
        self.MSI_only=MSI_only
        #MSI_only=True if want to only include samples with MSI determined
        self.PathologyList=PathologyList
        #PathologyList is a list of the keys from pathology dictionaries to be included in tileplot
        self.HRD_only=HRD_only
        self.exclude_VUS=exclude_VUS
        self.special_type_list=special_type_list


    def getVariantsInOncoPrintFormat(self):
        """
        Loads all the variants from the oneVariantPerLine.txt file.
        :param variants_path: A string with the path to the variants file
        """
        _log = logging.getLogger(name = 'LoadVariants')
        _log.info('Loading variants...')


        variantDict=self.varDict
        samples=self.specimens

        tile_variants=dict()
        line_number = 1
        status_interval = 100000
        if self.GroupDict:
            for group in self.GroupDict.keys():
                if group in self.GroupsToSkip: continue
                #for TRFid in self.GroupDict[group]["TRFlist"]:
                for TRFid in self.GroupDict[group]:
                    tile_variants[TRFid+"_"+ group]=[TRFid, "n/a", group, 11, "10000", "gene", "0", "0", "0"]

        for s in samples:
            line_number += 1
            if line_number%status_interval==0:
                _log.info('Processed ' + str(line_number) + ' variants...')

            MSI=s.msi_status
            trf=s.specimen_name
            disease=s.disease_ontology
            TMB=s.mutational_load_per_mb
            if self.HRD_only==True:
                HRD=s.HRD
                if HRD=="NA": continue

            #to print TMB and MSI status to include in tileplot
            if self.MSI_only==True:
                if MSI=="None" or MSI=="MSI unknown": continue
            if TMB == None:
                TMB = 9999999999999 # this number is > the length of the baitset, so is an impossible value to naturally occur
                TMBLevel="0"
            TMB=float(TMB)
            if TMB >=20 and TMB < 9999999999999: # this number is > the length of the baitset, so is an impossible value to naturally occur
                TMBLevel="3"
            elif TMB >=6 and TMB <20:
                TMBLevel="2"
            elif TMB <6:
                TMBLevel="1"
            if trf+"_"+"TMB" not in tile_variants :
                tile_variants[trf+"_"+"TMB"]=[trf, disease, "none", 10, TMB, "TMB", TMBLevel, "0", "0"]

            if MSI=="MSI-H":
                MSIstat="3"
            elif MSI=="MSI ambiguous":
                MSIstat="2"
            elif MSI=="MSS":
                MSIstat="1"
            else:
                MSIstat="4"
            if trf+"_"+"MSI" not in tile_variants:
                tile_variants[trf+"_"+"MSI"]=[trf, disease, "none", 0, "10000", "MSI", "0", MSIstat, "0"]

            #Scrap pathology data. Dexter recommends restricting it to internal FMI path staining for the sake of consistency in scoring. 
            if self.PathologyList:
                for PathType in self.PathologyList:
                    if PathType in s.pathology:
                        p=s.pathology[PathType]
                        result=p.status
                        testType=p.marker
                        if testType=="ER":
                            if result=="POSITIVE":
                                status="3"
                            elif result=="EQUIVOCAL":
                                status="2"
                            elif result=="NEGATIVE":
                                status="1"
                        elif testType=="PD-L1" or testType=="PD-1":
                            if result=="HIGH POSITIVE" or result=="MODERATE POSITIVE":
                                status="3"
                            elif result=="LOW POSITIVE":
                                status="2"
                            elif result=="NEGATIVE":
                                status="1"
                    else:
                        status="0"
                        testType = "PD-L1" #edit this to the specific pathology test of interest
                    if trf+"_"+testType not in tile_variants:
                        tile_variants[trf+"_"+testType]=[trf, disease, "none", 0, "10000", testType, "0", "0", status]

            if trf not in variantDict: continue
            variant=variantDict[trf]
            for v in variant:
                if self.exclude_VUS==True and v.driver_status_consensus==DriverStatus.unknown: continue
                variant_type=v.variant_type
                geneName=v.gene
                if variant_type=="RE":
                    geneName2=v.gene2
                    if geneName2 in self.GroupGeneList: continue
                if geneName in self.GroupGeneList: continue

                gene=''
                tile_type=''

                if variant_type=='SV':
                    (gene,tile_type) = self.parseSV(v)
                    if self.special_type_list:
                        if geneName + ":" + v.protein_effect in self.special_type_list:
                            tile_type= 11
                elif variant_type=='CN':
                    (gene,tile_type) = self.parseCN(v)
                elif variant_type=='RE':
                    (gene,tile_type) = self.parseRE(v)


                if '' not in (gene, tile_type):
                    if trf+"_"+gene in tile_variants and tile_variants[trf+"_"+gene][3]!=tile_type:
                        tile_types = sorted([tile_variants[trf+"_"+gene][3], tile_type])
                        if tile_types == [1,2]:  # all standard SV (not trunc, fs, or splice) + amp
                            tile_type = 6
                        elif tile_types == [3,4]:  # deleterious SVs + deletion + RE truncations
                            tile_type = 7
                        elif tile_types == [2,5]:  # amp + non-trunc RE
                            tile_type = 8
                        elif 11 in tile_types:
                            tile_type = 11
                        else:   # other multiple
                            tile_type = 9
                    tile_variants[trf+"_"+gene] = [trf, disease, gene, tile_type, "10000", "gene", "0", "0", "0"]
                    self.mutSamples.add(trf+'\t'+disease)
                elif trf not in self.mutSamples:
                    self.noMutSamples.add(trf+'\t'+disease) # if a sample is encountered for the first time with a mutation in a gene not in gene list, it will be included here.  Need to remove these later.




        self.noMutSamples = [sd for sd in self.noMutSamples if sd not in self.mutSamples]  # final check removing samples that are in self.mutSamples
        self.tile_variants=tile_variants

    def parseSV(self, v):
        gene = v.gene
        coding_type = v.coding_type
        if gene not in self.queryGenes:
            return ('', '')
        else:
            specialTypes = ['nonsense', 'frameshift','splice']
            if coding_type in specialTypes:
                tile_type = 4
            else:
                tile_type = 1

            return (gene,tile_type)

    def parseCN(self, v):
        gene = v.gene
        cn_alteration = str(v.amp_or_del)
        cn_alt=cn_alteration.split(".")[1]
        if gene not in self.queryGenes:
            return ('', '')
        else:
            toTileType_dict = {'amplification': 2, 'deletion': 3}
            tile_type = toTileType_dict[cn_alt]
            return (gene,tile_type)

    def parseRE(self, v):
        varGene = v.gene
        otherGene = v.gene2
        if varGene == 'KMT2A':
            varGene = 'MLL'
        type_re = v.coding_type


        if type_re == 'truncation':
            tile_type = 4
        else:
            tile_type = 5
        if varGene not in self.queryGenes:
            if otherGene not in self.queryGenes:
                return ('', '')
            else:
                return(otherGene, tile_type)
        return (varGene,tile_type)

    def outputOncoPrintVariants(self, outFile, geneoutfile):
        self.getVariantsInOncoPrintFormat()
        f=open(outFile,'w')
        header=['trf', 'disease', 'gene', 'type', 'TMB', 'functional', "TMBLevel", "MSILevel", "PathLevel"]
        f.write("\t".join(header)+'\n')

        # Convert the dict of samples/gene (sg) to a list of variant info lists to print
        variantsToPrint=list()
        for (sg,info) in self.tile_variants.iteritems():
            variantsToPrint.append(info)

        # a counter to track how often each gene appears
        counts = collections.Counter(t[2] for t in variantsToPrint)
        TotalSampleCount = collections.Counter(k[5] for k in variantsToPrint)

        SampleNum=float(TotalSampleCount["MSI"])

        #to give only top ten frequent genes
        #this is necessary for when querying all genes
        if self.TopTen==True:
            TopTenGenes={gene: counts[gene] for gene in sorted(counts, key=counts.get, reverse=True)[:11]}
            for key in counts.keys():
                value=counts[key]
                if key not in TopTenGenes.keys():
                    del counts[key]
        #this above should make it so only give tileplot genes based on an input list

        #to write out a geneList file
        g=open(geneoutfile, 'w')
        if self.TopTen==True:
            for gene1 in sorted(TopTenGenes, key=TopTenGenes.get, reverse=True):
                g.write(gene1+'\n')
            g.close()
        elif self.TopTen==False:
            for gene in counts.keys():
                g.write(gene+'\n')
            g.close()


        # order the list by gene frequency, then by tile_type
        for variantInfo in sorted(variantsToPrint, key=lambda x: (counts[x[2]], x[3], x[4], x[7]), reverse=True):
            if variantInfo[2] in counts.keys() or variantInfo[2]=="none":
                f.write("\t".join(map(str,variantInfo))+'\n')


        # Print samples with no mutations. Needed for correct denominator.
        for sd in self.noMutSamples:
            f.write(sd+"\tnone\t0\t10000\tgene\t0\t0\t0\n")
        f.close()
if __name__=="__main__":
    run(args)
