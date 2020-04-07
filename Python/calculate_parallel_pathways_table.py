from __future__ import division
import os, sys
import numpy as np

import phylo_tools as pt


MCR = 0.75

treatments = ['0','1','2']
taxa = ['B', 'C', 'D']
maple_types = ['signature', 'complex', 'pathway', 'function']


kegg_dict_count = {}
maple_dict_count = {}
maple_annotation_dict = {}

for treatment in treatments:

    #kegg_dict_count[treatment] = {}
    #maple_dict_count[treatment] = {}


    for taxon in taxa:

        protein_id_kegg_dict = {}

        protein_id_kegg = open(pt.get_path() + '/data/reference_assemblies_task2/MAPLE/%s_MAPLE_result/query.fst.ko' % taxon, 'r')
        # make protein ID => KEGG map
        for i, line in enumerate(protein_id_kegg):
            line = line.strip()
            items = line.split("\t")
            protein_id = items[0]
            if items[1] != 'K_NA':
                protein_id_kegg_dict[items[0]] = items[1]


        significant_genes = open(pt.get_path() + '/data/timecourse_final/parallel_genes_%s.txt' % (treatment+taxon) , 'r')

        first_line = significant_genes.readline()
        first_line = first_line.strip()
        first_line_items = first_line.split(",")

        kegg_list = []
        # get KEGG annotation of genes with significant multiplicity
        for i, line in enumerate(significant_genes):
            line = line.strip()
            items = line.split(",")
            protein_id = items[1].strip()
            if protein_id in protein_id_kegg_dict:
                kegg_list.append(protein_id_kegg_dict[protein_id])

                if protein_id_kegg_dict[protein_id] not in kegg_dict_count:
                    kegg_dict_count[ protein_id_kegg_dict[protein_id]] = {}

                if treatment not in kegg_dict_count[ protein_id_kegg_dict[protein_id]]:
                    kegg_dict_count[ protein_id_kegg_dict[protein_id]][treatment] = []

                if taxon not in kegg_dict_count[protein_id_kegg_dict[protein_id]][treatment]:
                    kegg_dict_count[protein_id_kegg_dict[protein_id]][treatment].append(taxon)


                #if taxon not in kegg_dict_count[treatment][protein_id_kegg_dict[protein_id]]:
                #    kegg_dict_count[treatment][protein_id_kegg_dict[protein_id]].append(taxon)

        # map KEGG genes onto pathays
        # get pathways
        # go through signatures, complexes, pathways, and functions
        maple_to_keep = {}
        for maple_type in maple_types:
            maple_file = open(pt.get_path() + '/data/reference_assemblies_task2/MAPLE/%s_MAPLE_result/module_%s.tsv' % (taxon, maple_type) , 'r')

            first_line_maple = maple_file.readline()
            first_line_maple = first_line_maple.strip()
            first_line_items_maple = first_line_maple.split("\t")

            for i, line in enumerate(maple_file):
                line = line.strip()
                items = line.split("\t")

                #if float(items[10]) < MCR:
                #    continue
                # remove rows that are less than 80% complete
                # query(coverage) = MCR % (ITR)
                # query(coverage/max) = MCR % (WC)
                # query(coverage/mode) = Q-value

                maple_name = items[0].split('_')[0]
                maple_to_keep[maple_name] = {}#[items[1], items[2]]
                maple_to_keep[maple_name]['type'] = items[1]
                maple_to_keep[maple_name]['description'] = items[2]

                if maple_name not in maple_annotation_dict:
                    maple_annotation_dict[maple_name] = {}
                    maple_annotation_dict[maple_name]['type'] = items[1]
                    maple_annotation_dict[maple_name]['description'] = items[2]

                maple_matrix = open(pt.get_path() + '/data/reference_assemblies_task2/MAPLE/%s_MAPLE_result/KAAS/%s_matrix.txt' % (taxon, maple_name) , 'r')
                maple_kegg = []
                # only get kegg genes that exist in the taxon's genome
                for matrix_line_i, matrix_line in enumerate(maple_matrix):
                    matrix_line = matrix_line.strip()
                    matrix_items = matrix_line.split('\t')
                    if 'K' not in matrix_items[0]:
                        continue
                    if len(matrix_items) < 3:
                        continue
                    if int(matrix_items[1]) == 0:
                        continue
                    maple_kegg.append(matrix_items[0])

                maple_to_keep[maple_name]['KEGG'] = maple_kegg

        #significant_kegg_maple_map = {}
        # now map kegg to pathways
        for kegg_i in kegg_list:

            for maple_i, maple_i_dict in maple_to_keep.items():

                if kegg_i in maple_i_dict['KEGG']:

                    #if maple_i not in maple_dict_count[treatment]:
                    #    maple_dict_count[treatment][maple_i] = []

                    #if taxon not in maple_dict_count[treatment][maple_i]:
                    #    maple_dict_count[treatment][maple_dict_count].append(taxon)


                    if maple_i not in maple_dict_count:
                        maple_dict_count[ maple_i] = {}

                    if treatment not in maple_dict_count[maple_i]:
                        maple_dict_count[ maple_i][treatment] = []

                    if taxon not in maple_dict_count[maple_i][treatment]:
                        maple_dict_count[maple_i][treatment].append(taxon)


                    #if kegg_i not in significant_kegg_maple_map:
                    #    significant_kegg_maple_map[kegg_i] = []

                    #significant_kegg_maple_map[kegg_i].append( maple_i)



convergence_table = open(pt.get_path() + '/data/parallel_pathways.txt', 'w')

# alot of ribosome complex, merge those first
convergence_table.write(", ".join(["Module ID", "Module type", "Description", "1-day", "10-day", "100-day" ]) + '\n' )

for maple_i, maple_i_dict in maple_dict_count.items():
    keep_maple=False
    for treatment in treatments:
        if treatment in  maple_i_dict:
            if len(maple_i_dict[treatment]) >= 2:
                keep_maple=True

    if keep_maple == False:
        continue

    type = maple_annotation_dict[maple_i]['type']
    description = maple_annotation_dict[maple_i]['description']

    if 'Ribosome' in description:
        continue

    if '0' in maple_i_dict:
        out_0 = str(len(maple_i_dict['0'])) + '/' + str(len(taxa))
    else:
        out_0 = '0/' + str(len(taxa))

    if '1' in maple_i_dict:
        out_1 = str(len(maple_i_dict['1'])) + '/' + str(len(taxa))
    else:
        out_1 = '0/' + str(len(taxa))

    if '2' in maple_i_dict:
        out_2 = str(len(maple_i_dict['2'])) + '/' + str(len(taxa))
    else:
        out_2 = '0/' + str(len(taxa))


    convergence_table.write(", ".join([maple_i, type, description, out_0, out_1, out_2 ]) + '\n' )

convergence_table.close()
#for maple_i, maple_i_list in maple_dict_count.items():
#    if len(maple_i_list) > 1:
#        print(maple_i, maple_annotation_dict[maple_i]['type'], maple_annotation_dict[maple_i]['description'],  maple_i_list )
