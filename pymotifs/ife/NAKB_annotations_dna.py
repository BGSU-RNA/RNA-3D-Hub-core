import csv

from pymotifs import core
from pymotifs import models as mod
from sqlalchemy import desc

class Loader(core.SimpleLoader):

    merge_data = True
    allow_no_data = True


    def to_process(self, pdbs, **kwrags):
        return ['dna']
        table_list = []
        # with open("/usr/local/pipeline/hub-core/NAKB_annotations_dna.csv", "r") as annotations:
        #     nakb_csv = csv.reader(annotations, delimiter=',')
            
        #     for row in nakb_csv:
        #         if row[0] != 'id':
        #             NA_dict = {}
        #             NA_dict['pdb_id'] = row[1]
        #             NA_dict['property'] = 'NAKB_NA_annotation'
        #             NA_dict['value'] = row[2]
        #             table_list.append(NA_dict)
        #         if row[0] != 'pdbid':
        #             prot_dict = {}
        #             print('second if')
        #             prot_dict['pdb_id'] = row[1]
        #             prot_dict['property'] = 'NAKB_protein_annotation'
        #             prot_dict['value'] = row[3]
        #             table_list.append(prot_dict)           
        # print(table_list[0])
        # print(table_list[1])
        # return table_list        

    def get_nakb_annotations(self):
        # table_list = []
        table = []
        with open("https://www.nakb.org/node/solr/nakb/select?fl=id,pdbid,NAKBnaList,NAKBprotList&q=NAKBna:*%20OR%20NAKBprot:*&wt=csv&rows=200000", "r") as annotations:
            nakb_csv = csv.reader(annotations, delimiter=',')
            
            for row in nakb_csv:
                if row[0] != 'id':
                    # NA_dict = {}
                    # NA_dict['pdb_id'] = row[1]
                    # NA_dict['property'] = 'NAKB_NA_annotation'
                    # NA_dict['value'] = row[2]
                    # table_list.append(NA_dict)
                    table.append(mod.PdbPropertyValue(
                                    pdb_id = row[1],
                                    property = 'NAKB_NA_annotation',
                                    value = row[2]))
                if row[1] != 'pdbid':
                    # prot_dict = {}
                    # prot_dict['pdb_id'] = row[1]
                    # prot_dict['property'] = 'NAKB_protein_annotation'
                    # prot_dict['value'] = row[3]
                    # table_list.append(prot_dict)     
                    table.append(mod.PdbPropertyValue(
                                    pdb_id = row[1],
                                    property = 'NAKB_protein_annotation',
                                    value = row[3]))      
        # print(table_list[0])
        # print(table_list[1])
        # for row in table_list:

        return table   

    def has_data(self, args, **kwargs):
        return False

    def query(self, session, pdb):
        return session.query(mod.dna_annotations)


    def data(self, annotation_row, **kwargs):
        table_list = self.get_nakb_annotations()
        # nakb_row = []
        # nakb_row.append(mod.PdbPropertyValue(
        #                 pdb_id = annotation_row['pdb_id'],
        #                 property = annotation_row['property'],
        #                 value = annotation_row['value']))
        return table_list
        # return mod.DnaAnnotations(annotation_row)

