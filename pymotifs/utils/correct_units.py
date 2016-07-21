"""This package contains tools to try and correct poorly translated unit ids.
"""

import collections as coll

import fr3d.unit_ids as fuid

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils import row2dict

Unit = coll.namedtuple('Unit', ['pdb', 'model', 'chain', 'component_number',
                                'insertion_code', 'alt_id', 'symmetry'])


class Correcter(core.Base):
    """This is a class to try and correct some of the mistakes that are present
    in our mappings from old to new style ids. Basically, this compares a unit
    id to the known unit ids. If there is not the same id some corrections are
    attempted to make the units match.
    """

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
                                  mod.UnitInfo.pdb_id.label('pdb'),
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
        del parts['atom_name']
        return Unit(**parts)

    def correct_nothing(self, unit):
        """Do nothing to unit.

        :param: The unit
        :returns: The same unit.
        """
        return Unit(*unit)

    def correct_model(self, unit):
        """If the model is not 1 and the symmetry operator is not 1_555 then
        create a unit with model as 1. If that is not true this will return
        None.

        :param: The unit
        """

        if unit.model != 1 and unit.symmetry != '1_555':
            return unit._replace(model=1)
        return None

    def correct_alt_id(self, unit):
        """If the alt id is None, then try create a new unit with with alt id
        'A', otherwise return None.
        """

        if unit.alt_id is None:
            return unit._replace(alt_id='A')
        return None

    def correct_model_and_alt_id(self, unit):
        """Alter a model and alt id if needed, if not returns None.
        """

        if unit.alt_id is None and unit.model != 1 and \
                unit.symmetry != '1_555':
            return unit._replace(alt_id = 'A', model=1)
        return None

    def correct(self, mapping, unit):
        """Attempt to correct a single unit.
        """

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

    def correct_structure(self, pdb, mapping, units, require_all=True):
        """Correct units from a single structure.

        :param str pdb: The PDB id to use.
        :param dict mapping: The mapping to use.
        :param list units: The list of Unit's to correct.
        :returns: A list of the corrected units.
        """

        valid = []
        for unit in units:
            if unit.pdb != pdb:
                continue

            fixed = self.correct(mapping, unit)
            if not fixed:
                msg = "Could not correct {unit}".format(unit=str(unit))
                if require_all:
                    raise core.InvalidState(msg)
                else:
                    self.logger.error(msg)
                continue
            valid.append(fixed)

        if len(valid) != len(set(valid)):
            raise core.InvalidState("Did not produce unique normalization")

        return valid

    def __call__(self, unit_ids, require_all=True):
        units = [self.as_unit(uid) for uid in unit_ids]
        pdbs = set(unit.pdb for unit in units)

        valid = []
        for pdb in pdbs:
            normalization = self.normalization(pdb)
            valid.extend(self.correct_structure(pdb, normalization, units,
                                                require_all=require_all))

        return valid


class TranslateCorrectly(core.Base):
    """Translate and correct some unit ids. This will use the translate tools
    to turn an old style to a new style id followed by attempting to correct
    the units id.
    """

    def correct(self, pdb, unit_ids):
        """Correct the units in a given pdb file.

        :param str pdb: The PDB id to use.
        :param list unit_ids: The unit ids to correct.
        :returns: The corrected unit ids.
        """

        correcter = Correcter(self.config, self.session)
        mapping = correcter.normalized_mapping(pdb)
        units = [correcter.as_unit(uid) for uid in unit_ids]
        corrected = correcter.correct_structure(pdb, mapping, units)
        return [mapping[u] for u in corrected]

    def __call__(self, pdb, unit_ids, sort=False):
        """Translate and correct unit ids in a given structure.

        :param str pdb: The PDB to use.
        :param str unit_ids: A comma seperated string of unit ids.
        :param bool sort: Weather or not we should sort the unit ids.
        :returns: A list sorted or not of the translated and corrected unit ids
        given.
        """

        translator = Translator(self.session)
        translated = translator(unit_ids)
        corrected = self.correct(pdb, translated)
        if sort:
            return sorted(corrected)
        return corrected
