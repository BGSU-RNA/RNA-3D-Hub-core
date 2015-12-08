"""This is the main module for loading integrated functional elements into the
database. Integrated functional elements (ife) are groups of 1 or more chains
that are used to cluster into the nr set. Each element is used as the basis for
nr clustering and analysis.

Chains within an element are in two different categories, integral or
accompanying. There is always one or more integral chain. The integral chain(s)
are structured chains which interact with themselves to form the functional
unit. They are permanently a part of the structure and required for its
function. For example, this could be a tRNA or the 5.8S and LSU. In the second
case the LSU and 5.8S function as one unit and we attempt to recognize and
respect with that in our nr groupings.

Accompanying chains however, are chains that are not permanent and not required
for the function. They may also be accompanying because they lack a 3D
structure. For example, an mRNA that has no internal structure and that
particular mRNA is not required for tRNA function.

Detection of integral vs accompanying is done on the basis of the extent of the
interactions between the two chains. For details on that logic look at
pymotifs.ife.grouper.
"""

from pymotifs import core

from pymotifs.ife.info import Loader as InfoLoader


class Loader(core.StageContainer):
    """This loader will load all ife related data into the database. It is a
    mass loader so it will load the information and then load the chain
    specific information. For details on what each step does look at
    pymotifs.ife.info and pymotifs.ife.chains.
    """

    stages = set([InfoLoader])
