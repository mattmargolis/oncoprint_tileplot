#!/usr/bin/env python

from researchmaster.report import *
import argparse

import getDataForTilePlot_Modular
reload(getDataForTilePlot_Modular) 
from getDataForTilePlot_Modular import OncoPrinter

###############################################################################
# Main command line script for performing "rm pulls", i.e. queries of the
# Research Master freezes into excel documents for delivery to customers.
#
# Author: Tim Fennell
###############################################################################

# Configure logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s %(levelname)s %(name)s] %(message)s', datefmt='%m/%d/%Y %H:%M:%S')

def enum_converter(enum):
    """Function that builds and returns a function that converts from String to a specific enum type."""
    def convert(value):
        try:
            return enum[value]
        except:
            raise argparse.ArgumentTypeError('{0} is not a valid {1}'.format(value, enum))
    return convert

###############################################################################
# Arg parsing
###############################################################################
e_to_n = lambda enum_value: enum_value.name
dummy = QueryParams()

parser = argparse.ArgumentParser(description='Run a query against research master and produce excel reports.',
                                 fromfile_prefix_chars='@')
# non-query params
parser.add_argument('--freeze', '-f', metavar='freeze_dir', help='path to the freeze directory to query', required=True)
parser.add_argument('--aux',    '-x', metavar='aux_dir', help='aux directory, defaults to freeze directory', required=False)
parser.add_argument('--out',    '-o', metavar='file', help='base name for output files', required=True)
# query params
parser.add_argument('--genes',    metavar='gene', nargs='*', default=[], help='genes to query for all types of variant')
parser.add_argument('--sv-genes', metavar='gene', nargs='*', default=[], help='genes to query for short variants')
parser.add_argument('--cn-genes', metavar='gene', nargs='*', default=[], help='genes to query for copy number variants')
parser.add_argument('--re-genes', metavar='gene', nargs='*', default=[], help='genes to query for rearrangements')
parser.add_argument('--nh-genes', metavar='gene', nargs='*', default=[], help='genes to query for non human variants')
parser.add_argument('--negative-genes', metavar='gene', nargs='*', default=[],
                    help='eliminate samples with known/likely variants in genes')
parser.add_argument('--bait-sets', metavar='baitset', nargs='*', type=enum_converter(BaitSet),
                    default=dummy.bait_sets, help='One or more of {0}. Default: all.'.format(map(e_to_n, BaitSet.all())))
parser.add_argument('--driver-statuses', metavar='status', nargs='*', type=enum_converter(DriverStatus),
                    default=dummy.driver_statuses,
                    help='One or more of {0}. Default: {1}.'.format(map(e_to_n, DriverStatus.all()), map(e_to_n, dummy.driver_statuses)))
parser.add_argument('--diseases', metavar='disease', nargs='*', default=[], help='diseases to query for')
parser.add_argument('--negative-diseases', metavar='disease', nargs='*', default=[], help='diseases to exclude')
parser.add_argument('--min-samples-per-disease', type=int, default=dummy.min_samples_per_disease,
                    help='exclude diseases from query if they have fewer than this many assayed samples. Default: {0}'.format(dummy.min_samples_per_disease))
parser.add_argument('--primary-met-statuses', metavar='status', nargs='*', type=enum_converter(PrimaryMetStatus), default=dummy.primary_met_statuses,
                    help='One or more of {0}. Default {1}.'.format(map(e_to_n, PrimaryMetStatus.all()), map(e_to_n, dummy.primary_met_statuses)))
parser.add_argument('--trfs', nargs='*', default=[], help="specific sample TRFs to query for")
parser.add_argument('--xrns', nargs='*', default=[], help="specific sample XRNs to query for")
parser.add_argument('--min-age', type=int, help='minimum age of patient to restrict query to')
parser.add_argument('--max-age', type=int, help='maximum age of patient to restrict query to')
parser.add_argument('--studies', metavar='studies', nargs='*', default=[], help='Study filter')
parser.add_argument('--genes_sv_maf_min', type=float, help='minimum maf threshold for short variants')
parser.add_argument('--genes_sv_maf_max', type=float, help='maximum maf threshold for short variants')
parser.add_argument('--computational_tumor_purity_min', type=float, help='minimum computational tumor purity threshold')
parser.add_argument('--computational_tumor_purity_max', type=float, help='maximum computational tumor purity threshold')
# non-query params 
parser.add_argument('--report-queried-specimens', action="store_true", default=False,
                    help='if True, report writes queried_specimens tab to include samples without specified variants)')
parser.add_argument('--filter_obsolete_off', action="store_true", default=False,
                    help='if True, turns filter_obsolete off, allowing non-clinical studies and non-default (T5a, T7, D2, CF2) baitsets')
parser.add_argument('--filter_multiples_off', action="store_true", default=False,
                    help='if True, turns filter_multiples off, allowing multiple specimens from same patient')
parser.add_argument('--filter_nonresearch_off', action="store_true", default=False,
                    help='if True, turns filters for research permissions off (obtain necessary approval separately)')
parser.add_argument('--filter_duplicates', action="store_true", default=False,
                    help='if True, turns filter_duplicates on, preventing any sample with the same TRF listed multiple times from being included')
parser.add_argument('--skip_loading_variants', action="store_true", default=False,
                    help='if True, skips loading variants into freeze')
parser.add_argument('--skip_loading_vus', action="store_true", default=False,
                    help='if True, skips loading VUS into freeze')

###############################################################################
# Parse the args and run the query
###############################################################################
args = parser.parse_args()

if not args.aux:
    args.aux = args.freeze

freeze_params=FreezeParams(bait_sets=args.bait_sets, 
                           diseases=args.diseases,
                           skip_loading_variants=args.skip_loading_variants,
                           skip_loading_vus=args.skip_loading_vus)

params = QueryParams(genes=args.genes,
                     genes_sv=args.sv_genes,
                     genes_cn=args.cn_genes,
                     genes_re=args.re_genes,
                     genes_nh=args.nh_genes,
                     negative_genes=args.negative_genes,
                     bait_sets=args.bait_sets,
                     driver_statuses=args.driver_statuses,
                     diseases=args.diseases,
                     negative_diseases=args.negative_diseases,
                     min_samples_per_disease=args.min_samples_per_disease,
                     specimen_trfs=args.trfs,
                     specimen_xrns=args.xrns,
                     min_age=args.min_age,
                     max_age=args.max_age,
                     primary_met_statuses=args.primary_met_statuses,
                     studies=args.studies,
                     computational_tumor_purity_min=args.computational_tumor_purity_min,
                     computational_tumor_purity_max=args.computational_tumor_purity_max,
                     genes_sv_maf_min=args.genes_sv_maf_min,
                     genes_sv_maf_max=args.genes_sv_maf_max
                     )

logging.info('Loading freeze...')
if args.studies:
    logging.info('Study list was detected. Study samples may require the --filter_obsolete_off, --filter_multiples_off, and --filter_nonresearch_off flags to all be enacted to be queryable...')
if args.filter_obsolete_off:
    filter_obsolete=False
else:
    filter_obsolete=True
if args.filter_multiples_off:
    filter_multiples=False
else:
    filter_multiples=True
if args.filter_nonresearch_off:
    filter_exUS=False
else:
    filter_exUS=True
freeze = Freeze(freeze_path=args.freeze, aux_path=args.aux, filter_obsolete=filter_obsolete, filter_multiples=filter_multiples, filter_duplicates=args.filter_duplicates, filter_exUS=filter_exUS, freeze_params=freeze_params)

logging.info('Querying freeze...')
query = Query(freeze, params)
query.execute()

logging.info('Writing reports...')
for raw, fn in [(True, 'raw'), (False, 'semi-final')]:
    filename = args.out + '.' + fn + '.xlsx'
    writer = ExcelWriter(query=query, filename=filename, raw=raw, report_queried_specimens=args.report_queried_specimens)
    writer.write()

logging.info('Finished.')






# genesUsed = [list of genes wanted to make tileplot]
# driver_string I use 'KL' as a driver_string here 

genesUsed = query.params.genes_sv | query.params.genes_cn | query.params.genes_re | query.params.genes_nh
print genesUsed

# TILE PLOT VIS DATA GEN
specimen_name_to_variants = defaultdict( lambda: [] )
for s, vs in groupby( sorted( query.found_variants, key=lambda v: v.specimen_name ), key=lambda v: v.specimen_name ):
    specimen_name_to_variants[s] = list( vs )

specialDict = dict()

oncoprinter = OncoPrinter( varDict=specimen_name_to_variants, specimens=query.queried_specimens, queryGenes=genesUsed,
                          GroupDict=specialDict, GroupGeneList=[], GroupsToSkip=[], TopTen=False, MSI_only=False,
                          PathologyList=[], HRD_only=False, exclude_VUS=False, special_type_list=[] )

outfileTile = args.out + '.tileplotdata.txt'
geneoutfileTile = args.out +'.gene_list.txt'
oncoprinter.outputOncoPrintVariants( outfileTile, geneoutfileTile ) 




 # COMMENTED OUT
 # out_tile_file = './vis/%s_%s.raw.xlsx' % ( out_base, driver_string )

 # ( out_path, excel_out_file ) = os.path.split( out_tile_file ) # run this to make sure

 # remove both the .xlsx then the .raw extensions from the filename
 # basefile = os.path.splitext( os.path.splitext( excel_out_file )[0])[0]
 # outfileTile = os.path.join( out_path,basefile + '.tileplotdata.txt' )
 # geneoutfileTile = os.path.join( out_path,basefile+'.gene_list.txt' )
 # oncoprinter.outputOncoPrintVariants( outfileTile, geneoutfileTile ) # this will write the files
 # print driver_string