"""This package contains tools to try and correct poorly translated unit ids.
"""

import collections as coll

import fr3d.unit_ids as fuid

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict

Unit = coll.namedtuple('Unit', ['model', 'chain', 'component_number',
                                'insertion_code', 'alt_id', 'symmetry'])


class Correcter(core.Base):
    def normalized_mapping(self, pdb_id):
        """This produces a dictonary that can be used to correct bad unit ids.
        Some of the loops stored after we migrated the database have incorrect
        unit ids. The errors appear to be of 2 kinds, incorrect model number
        and possibly bad alt ids. By producing this mapping we try to correct
        the issue by finding the correct unit id.

        :param str pdb_id: The pdb id to get units for.
        :returns: A dictonary with Unit keys mapping to the unit id.
        """
        with self.session() as session:
            query = session.query(mod.UnitInfo.unit_id,
                                  mod.UnitInfo.model,
                                  mod.UnitInfo.chain,
                                  mod.UnitInfo.number.label('component_number'),
                                  mod.UnitInfo.ins_code.label('insertion_code'),
                                  mod.UnitInfo.alt_id,
                                  mod.UnitInfo.sym_op.label('symmetry'),
                                  ).\
                filter(mod.UnitInfo.pdb_id == pdb_id)

            if not query.count():
                raise core.InvalidState("No units in %s" % pdb_id)

            mapping = {}
            for result in query:
                data = row2dict(result)
                unit_id = data.pop('unit_id')
                key = Unit(**data)
                if key in mapping:
                    raise core.InvalidState("Non unique mapping")
                mapping[key] = unit_id
        return mapping

    def as_unit(self, unit_id):
        """Turn a unit id into a Unit.

        :param str unit_id: The unit id.
        """

        parts = fuid.decode(unit_id)
        del parts['component_id']
        del parts['pdb']
        del parts['atom_name']
        return Unit(**parts)

    def correct_nothing(self, unit):
        return unit

    def correct_model(self, unit):
        if unit.model != 1 and unit.symmetry != '1_555':
            return unit._replace(model=1)
        return None

    def correct_alt_id(self, unit):
        if unit.alt_id is None:
            return unit._replace(alt_id='a')
        return None

    def correct_model_and_alt_id(self, unit):
        if unit.alt_id is None and unit.model != 1 and \
                unit.symmetry != '1_555':
            return unit._replace(alt_id = 'a', model=1)
        return None

    def correct(self, mapping, unit):
        corrections = [self.correct_nothing, self.correct_model,
                       self.correct_alt_id, self.correct_model_and_alt_id]
        for correction in corrections:
            norm = correction(unit)
            if norm:
                self.logger.debug("Attempting %s", correction.__name__)
                if norm in mapping:
                    self.logger.debug("Attempt worked")
                    return norm
                self.logger.debug("Attempt failed")
        return None

    def __call__(self, unit_ids, require_all=True):
        units = [self.as_unit(uid) for uid in unit_ids]
        pdbs = set(unit.pdb for unit in units)

        valid = []
        for pdb in pdbs:
            normalization = self.normalization(pdb)
            for unit in units:
                if unit.pdb != pdb:
                    continue

                fixed = self.correct(unit)
                if not fixed:
                    msg = "Could not correct %s" % unit
                    if require_all:
                        raise core.InvalidState(msg)
                    else:
                        self.logger.error(msg)
                valid.append(fixed)

        if len(valid) != len(set(valid)):
            raise core.InvalidState("Did not produce unique normalization")

        return valid
