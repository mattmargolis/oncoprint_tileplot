import getDataForTilePlot_Modular
reload(getDataForTilePlot_Modular) 
from getDataForTilePlot_Modular import OncoPrinter


# genesUsed = [list of genes wanted to make tileplot]
# driver_string I use 'KL' as a driver_string here 

# TILE PLOT VIS DATA GEN
    specimen_name_to_variants = defaultdict( lambda: [] )
    for s, vs in groupby( sorted( query.found_variants, key=lambda v: v.specimen_name ), key=lambda v: v.specimen_name ):
        specimen_name_to_variants[s] = list( vs )

    specialDict = dict()

    oncoprinter = OncoPrinter( varDict=specimen_name_to_variants, specimens=query.queried_specimens, queryGenes=genesUsed,
                          GroupDict=specialDict, GroupGeneList=[], GroupsToSkip=[], TopTen=False, MSI_only=False,
                          PathologyList=[], HRD_only=False, exclude_VUS=False, special_type_list=[] )

    out_tile_file = './vis/%s_%s.raw.xlsx' % ( out_base, driver_string )

    ( out_path, excel_out_file ) = os.path.split( out_tile_file ) # run this to make sure

    # remove both the .xlsx then the .raw extensions from the filename
    basefile = os.path.splitext( os.path.splitext( excel_out_file )[0])[0]
    outfileTile = os.path.join( out_path,basefile + '.tileplotdata.txt' )
    geneoutfileTile = os.path.join( out_path,basefile+'.gene_list.txt' )
    oncoprinter.outputOncoPrintVariants( outfileTile, geneoutfileTile ) # this will write the files
    print driver_string
