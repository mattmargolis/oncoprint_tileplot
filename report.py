from collections import OrderedDict
from query import *
from os.path import basename
import sys
import inspect
import os

version = 'Q4_2018_Freeze'

"""
Module that holds classes and functions for writing an RmQuery object out into excel files
"""

import xlsxwriter

class PlusPlus(object):
    """
    Stupid little counter class because there is no i++
    """
    def __init__(self, initial_value=0):
        self.i = initial_value

    def nxt(self):
        tmp = self.i
        self.i += 1
        return tmp

    def get(self):
        return self.i

    def reset(self, to_value=0):
        self.i = to_value


class ExcelWriter(object):
    __version__ = version

    def __init__(self, query, filename, raw=False, report_queried_specimens=False):
        self.query = query
        self.filename = filename
        self.raw = raw
        self.report_queried_specimens = report_queried_specimens
        
        self.workbook = xlsxwriter.Workbook(filename,{'nan_inf_to_errors': True})
        self.header_style  = self.workbook.add_format({'bold': True})
        self.wrapped_style = self.workbook.add_format({'text_wrap' : True, 'align': 'left'})
        self.right_align_style = self.workbook.add_format({'align': 'right'})
        self.percent_style = self.workbook.add_format({'num_format': '0.00%'})
        
    def write(self, extra_var_col = None, specify_var_col = None, extra_sample_col = None, specify_sample_col = None):
        if self.raw: self._write_manifest_tab(self.workbook)
        self._write_summary_tab(self.workbook, self.raw)
        self._write_variants_tab(self.workbook, self.raw, extra_var_col, specify_var_col)
        self._write_samples_tab(self.workbook, self.raw, extra_sample_col, specify_sample_col)
        self._write_diseases_tab(self.workbook, self.raw)
        self._write_groups_tab(self.workbook, self.raw)
        self._write_extra_tabs(self.workbook, self.raw)
        self.workbook.close()


    def list_string(self, list, sep=', '):
        """Helper function to turn a list of things into a comma separated string, and limit length."""
        string = sep.join(list)
        if len(string) > 32767:
            string = string[0:32764] + "..."
        return string

    def hugo_translation_summary(self):
        """ Helper function to write out """
        query = self.query
        ps = query.params
        queryName_translatedName = []
        if ps.is_query_all_genes():
            return "All"
        else:
            for x in query.original_query_genes:
                if x != query.freeze.hugo_gene_lookup[x]:
                    queryName_translatedName.append(x + ":" + query.freeze.hugo_gene_lookup[x])
            return(queryName_translatedName)


    def _write_manifest_tab(self, wb):
        """
        Writes a 'manifest' tab that records various information about the freeze used for this query/report
        and the command line and parameters used, so that it can easily be reconstructed later.
        """
        query = self.query
        freeze = query.freeze

        sheet = wb.add_worksheet("Manifest")
        row = PlusPlus(initial_value=1)
        col = PlusPlus()

        # Write out the headers
        sheet.write(0, col.nxt(), "filename", self.header_style)
        sheet.write(0, col.nxt(), "line_count", self.header_style)
        sheet.write(0, col.nxt(), "md5", self.header_style)

        # Write out the data
        for ff in freeze.freeze_files:
            rownum = row.nxt()
            col.reset()
            sheet.write(rownum, col.nxt(), ff.name)
            sheet.write(rownum, col.nxt(), ff.line_count)
            sheet.write(rownum, col.nxt(), ff.md5)

        ###Prints modified date for each Freeze/Query/Report file to help determine version used
        facts = OrderedDict()
        facts['Freeze Modification Date:'] = datetime.datetime.fromtimestamp(os.path.getmtime(inspect.getfile(Freeze))).strftime('%Y-%m-%d %H:%M:%S')
        facts['Query Modification Date:'] = datetime.datetime.fromtimestamp(os.path.getmtime(inspect.getfile(Query))).strftime('%Y-%m-%d %H:%M:%S')
        facts['Report Modification Date:'] = datetime.datetime.fromtimestamp(os.path.getmtime(inspect.getfile(ExcelWriter))).strftime('%Y-%m-%d %H:%M:%S')

        try:
            facts['Freeze Version:'] = Freeze.__version__
        except:
            facts['Freeze Version:'] = 'None'
        try:
            facts['Query Version:'] = Query.__version__
        except:
            facts['Query Version:'] = 'None'
        try:
            facts['Report Version:'] = ExcelWriter.__version__
        except:
            facts['Report Version:'] = 'None'

        for k, v in facts.items():
            rownum = row.nxt()
            sheet.write(rownum, 0, k, self.header_style)
            sheet.write(rownum, 1, v)


    def _write_summary_tab(self, wb, raw):
        """
        :type wb: xlsxwriter.Workbook
        """
        query = self.query
        freeze = query.freeze
        ps = query.params
        
        if ps.bait_sets == set(BaitSet.all()) and freeze.filter_obsolete == True:
            for obsolete_bait_set in freeze.bait_sets_to_filter:
                ps.bait_sets.remove(obsolete_bait_set)

        facts = OrderedDict()
        facts['Freeze Name']       = basename(freeze.freeze_path)
        facts['Bait Sets']                = self.list_string(map(lambda b: b.name,  ps.bait_sets))
        facts['Total Samples Queried']    = len(query.queried_specimens)
        facts['Variant Samples Reported'] = len(query.found_specimens)
        facts['Total Variants Reported']  = len(query.found_variants)
        facts['Unique Variants Reported'] = len(set(map(lambda v: v.long_id(), query.found_variants)))
        facts['Disease Ontologies']       = self.list_string(sorted(ps.diseases)) if ps.diseases else "All"
        facts['Genes Queried for Short Variants']       = "All" if ps.is_query_all_genes() else self.list_string(ps.genes_sv)
        facts['Genes Queried for Copy Number Variants'] = "All" if ps.is_query_all_genes() else self.list_string(ps.genes_cn)
        facts['Genes Queried for Rearrangements']       = "All" if ps.is_query_all_genes() else self.list_string(ps.genes_re)
        facts['Genes Queried for Non-Human Variants']       = "N/A" if ps.is_query_all_genes() else self.list_string(ps.genes_nh)
        if ps.negative_genes:
            facts['Excluding Samples Positive for Mutations In'] = self.list_string(ps.negative_genes)
        if raw:
            facts[' '] = ''
            facts['Freeze Params Used'] = 'False' if freeze.freeze_params is None else 'True'
            facts['Freeze Params: Baitsets'] = 'NA' if freeze.freeze_params is None else str(map(lambda x: x.name, freeze.freeze_params.bait_sets))
            facts['Freeze Params: Diseases'] = 'NA' if freeze.freeze_params is None else str(list(freeze.freeze_params.diseases))        
            facts['Freeze Params: Skipped Variants'] = 'NA' if freeze.freeze_params is None else str(freeze.freeze_params.skip_loading_variants)   
        if raw:
            facts['']             = ''
            facts['Command Line'] = " ".join(sys.argv)
            facts.update(ps.to_dict())

        facts['HUGO translated genes [GeneName:HUGO_GeneName]'] = "All" if self.hugo_translation_summary() == 'All' else self.list_string(self.hugo_translation_summary())
        sheet = wb.add_worksheet("Summary")
        row = 0
        for k, v in facts.items():
            if k != '' and k != ' ': k += ':'
            sheet.write(row, 0, k, self.header_style)
            sheet.write(row, 1, v, self.wrapped_style)
            row += 1

        sheet.set_column(0, 0, 33)
        sheet.set_column(1, 1, 120)

    def _write_variants_tab(self, wb, raw, extra_var_col = None, specify_var_col = None):
        """
        :type wb: xlsxwriter.Workbook
        :type extra_var_col: list; containing tuples of column name : attributes to be added to standard outputs
            if there are nested objects the attributes are listed in order, what was v.specimen.xrn.name is now ('header_name', 'specimen', 'xrn', 'name')
        :type specify_var_col: list; containing tuples of column name : attributes to replace standard outputs
        """
        var_col = [('xrn', 'specimen', 'xrn'),
                   ('trf', 'specimen_name'),
                   ('specimen_number', 'specimen_number'),
                   ('tissue', 'specimen', 'tissue'),
                   ('age [0-89]', 'specimen', 'patient_age'),
                   ('gender', 'specimen', 'gender', 'name'),
                   ('study', 'study_run_sheet'),
                   ('bait_set', 'bait_set', 'name'),
                   ('disease_ontology_term', 'disease_ontology'),
                   ('group_ontology_term', 'specimen', 'disease_group'),
                   ('gene', 'gene'),
                   ('partner_gene', 'gene2'),
                   ('pos1_RE', 'pos1'),
                   ('pos2_RE', 'pos2'),
                   ('comment_RE', 'comment'),
                   ('alteration_type', 'variant_type'),
                   ('cds_effect', 'transcript_effect'),
                   ('protein_effect', 'protein_effect'),
                   ('copy_number', 'copy_number'),
                   ('nonhuman_reads_per_million', 'RPM'),
                   ('fraction_reads', 'allele_freq'),
                   ('var_status', 'driver_status_consensus', 'name'),
                   ('long_id', 'long_id'),
                   ('Chromosome', 'chrom'),
                   ('position', 'position'),
                   ('refSeq_SV', 'ref_seq'),
                   ('altSeq_SV', 'alt_seq'),
                   ('transcript_SV', 'transcript'),
                   ('coding_type', 'coding_type')]
        if not raw:
            var_col = filter(lambda x: x[0] not in {'trf', 'specimen_number', 'study'}, var_col)
            
        if extra_var_col:
            var_col = var_col + extra_var_col
            
        if specify_var_col:
            var_col = specify_var_col
        
        query = self.query
        sheet = wb.add_worksheet("Variants")
        row = PlusPlus(initial_value=1)
        col = PlusPlus()

        # Write out the headers

        for column in var_col:
            sheet.write(0, col.nxt(), column[0], self.header_style)

        # Write out the data rows
        for v in query.found_variants:
            rownum = row.nxt()
            col.reset()

            for column in var_col:
                output = self._get_object_attribute(v, column[1:])
                try:
                    if column == ('comment_RE', 'comment') and output != '':
                        sheet.write(rownum, col.nxt(), output.decode('ascii'))
                    else:
                        sheet.write(rownum, col.nxt(), output)
                except UnicodeDecodeError:
                    sheet.write(rownum, col.get()-1, ''.join(i for i in output if ord(i) < 128))
                except TypeError:
                    sheet.write(rownum, col.get()-1, output())
            
        sheet.autofilter(0, 0, row.get()-1, col.get()-1)
        sheet.freeze_panes(1, 0)

    def _get_object_attribute(self, obj, attr_list):
        output = getattr(obj, attr_list[0], '')
        if len(attr_list) > 1:
            for attr in attr_list[1:]:
                output = getattr(output, attr, '')
        return output

    def _write_samples_tab(self, wb, raw, extra_sample_col = None, specify_sample_col = None):
        """:type wb: xlsxwriter.Workbook"""

        query = self.query

        sheet = wb.add_worksheet("Queried Samples" if self.report_queried_specimens else "Variant Samples")
        row = PlusPlus(initial_value=1)
        col = PlusPlus()

        sample_col = [('xrn', 'xrn'),
                      ('trf', 'specimen_name'),
                      ('specimen_number', 'specimen_number'),
                      ('disease_ontology_term', 'disease_ontology'),
                      ('group_ontology_term', 'disease_group'),
                      ('tissue', 'tissue'),
                      ('local_met_status', 'primary_or_met', 'name'),
                      ('age [0-89]', 'patient_age'),
                      ('gender', 'gender', 'name'),
                      ('pathologyPercentTumorNuclei', 'pathology_percent_tumor_nuclei'),
                      ('computationalTumorPurity', 'computational_tumor_purity'),
                      ('msi_status', 'msi_status'),
                      ('msi_pc1', 'msi_pc1'),
                      ('bait_set', 'bait_set', 'name'),
                      ('study', 'study_run_sheet'),
                      ('percent_genome_loh', 'loh_genome'),
                      ('mutation_load', 'mutational_load'),
                      ('mutation_load_per_mb', 'mutational_load_per_mb'),
                      ('short_variants', 'get_sv_summary'),
                      ('copy_number_alts', 'get_cn_summary'),
                      ('rearrangements', 'get_re_summary'),
                      ('non_human_variants', 'get_nh_summary')]
        
        if not raw:
            sample_col = filter(lambda x: x[0] not in {'trf', 'specimen_number', 'study'}, sample_col)
        
        if extra_sample_col:
            sample_col = sample_col + extra_sample_col
        
        if specify_sample_col:
            sample_col = specify_sample_col

        # Write out the headers

        for column in sample_col:
            sheet.write(0, col.nxt(), column[0], self.header_style)

        # Write out the values
        for s in (query.queried_specimens if self.report_queried_specimens else query.found_specimens):
            rownum = row.nxt()
            col.reset()

            for column in sample_col:
                output = self._get_object_attribute(s, column[1:])
                try:
                    sheet.write(rownum, col.nxt(), output)
                except TypeError:
                    sheet.write(rownum, col.get()-1, output())

        sheet.autofilter(0, 0, row.get()-1, col.get()-1)
        sheet.freeze_panes(1, 0)

    def _write_diseases_tab(self, wb, raw):
        """:type wb: xlsxwriter.Workbook"""

        sheet = wb.add_worksheet("Diseases")
        row = PlusPlus(initial_value=1)
        col = PlusPlus()
        
        query = self.query
        freeze = query.freeze
        ps = query.params
        
        if ps.bait_sets == set(BaitSet.all()) and freeze.filter_obsolete == True:
            for obsolete_bait_set in freeze.bait_sets_to_filter:
                ps.bait_sets.remove(obsolete_bait_set)

        # Write out the headers
        sheet.write(0, col.nxt(), "disease_ontology", self.header_style)
        sheet.write(0, col.nxt(), "group_ontology", self.header_style)
        sheet.write(0, col.nxt(), "variant_samples", self.header_style)
        sheet.write(0, col.nxt(), "total_samples", self.header_style)
        sheet.write(0, col.nxt(), "percent_variant", self.header_style)
        for bait_set in ps.bait_sets:
            sheet.write(0, col.nxt(), bait_set.name + " samples", self.header_style)

        # Write out the values
        for md in sorted(query.disease_metadata.values(), key=lambda md: md.disease_name):
            if md.disease_name in freeze.diseases_by_name:
                disease_group = freeze.diseases_by_name[md.disease_name].getGroup()
                group_name = ("underspecified" if disease_group is None else disease_group.name)
            else:
                group_name = "none"
            rownum = row.nxt()
            col.reset()
            sheet.write(rownum, col.nxt(), md.disease_name)
            sheet.write(rownum, col.nxt(), group_name)
            sheet.write(rownum, col.nxt(), md.found_count)
            sheet.write(rownum, col.nxt(), md.total())
            sheet.write(rownum, col.nxt(), md.found_count/float(md.total()), self.percent_style)
            for bait_set in ps.bait_sets:
                sheet.write(rownum, col.nxt(), md.found_by_bait_set[bait_set])

        sheet.autofilter(0, 0, row.get()-1, col.get()-1)
        sheet.freeze_panes(1, 0)

    def _write_groups_tab(self, wb, raw):
        """:type wb: xlsxwriter.Workbook"""

        sheet = wb.add_worksheet("Groups")
        row = PlusPlus(initial_value=1)
        col = PlusPlus()
        
        query = self.query
        freeze = query.freeze
        ps = query.params
        
        if ps.bait_sets == set(BaitSet.all()) and freeze.filter_obsolete == True:
            for obsolete_bait_set in freeze.bait_sets_to_filter:
                ps.bait_sets.remove(obsolete_bait_set)

        # Write out the headers
        sheet.write(0, col.nxt(), "disease_group", self.header_style)
        sheet.write(0, col.nxt(), "variant_samples", self.header_style)
        sheet.write(0, col.nxt(), "total_samples", self.header_style)
        sheet.write(0, col.nxt(), "percent_variant", self.header_style)
        for bait_set in ps.bait_sets:
            sheet.write(0, col.nxt(), bait_set.name + " samples", self.header_style)

        # Write out the values
        for md in sorted(query.group_metadata.values(), key=lambda md: md.disease_name):
            if md.disease_name in freeze.diseases_by_name and not freeze.diseases_by_name[md.disease_name].is_group:
                continue
            rownum = row.nxt()
            col.reset()
            sheet.write(rownum, col.nxt(), md.disease_name)
            sheet.write(rownum, col.nxt(), md.found_count)
            sheet.write(rownum, col.nxt(), md.total())
            sheet.write(rownum, col.nxt(), md.found_count/float(md.total()), self.percent_style)
            for bait_set in ps.bait_sets:
                sheet.write(rownum, col.nxt(), md.found_by_bait_set[bait_set])

        sheet.autofilter(0, 0, row.get()-1, col.get()-1)
        sheet.freeze_panes(1, 0)

    def _write_extra_tabs(self, wb, raw):
        """Does nothing, but allows sub-classes to override and create extra tabs."""
        None
