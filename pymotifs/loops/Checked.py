#!/usr/bin/env python
# coding: utf-8

# In[ ]:


""" Program to indicate whether loops have been checked or not. """

from pymotifs import core
from pymotifs import models as mod 

from pymotifs.loops.extractor import Loader as InfoLoader
from pymotifs.loops.positions import Loader as PositionLoader
from pymotifs.loops.quality import Loader as QualityLoader


class Loader(core.SimpleLoader):
    
    allow_no_data = True
    
    merge_data = True
    
    dependencies = set([InfoLoader, PositionLoader, QualityLoader])
    
    def table(self):
        return mod.PdbInfo
    
    def list_of_checked(self, pdbs, **kwargs):
        with self.session as session():
            query = session.query(mod.PdbInfo.pdb_id).                filter(mod.PdbInfo.loops_checked == 0)
            temp = [r.pdb_id for r in query]
            
        pdbs_checked = set(pdbs).intersection(temp)
            
        return pdbs_checked
        
    def data(self, pdbs, **kwargs):
        
        data = []
        
        for pdb in self.list_of_checked(pdbs):
            data.append({'pdb_id': pdb,
                        'loops_checked': 1
                        })
        return data

