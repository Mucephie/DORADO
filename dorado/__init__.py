# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *   # noqa
# ----------------------------------------------------------------------------

# __all__ = []
# from .example_mod import *   # noqa
from .fournax import * # noqa
from .vizual import * # noqa
from .rchive import * # noqa
from .waldo import * # noqa
from .imreduc import * # noqa
from .lazyeye import * # noqa

# __all__ += fournax.__all__
# __all__ += vizual.__all__
# __all__ += rchive.__all__
# __all__ += waldo.__all__
# __all__ += imreduc.__all__
# __all__ += lazyeye.__all__
# __all__ += dorado.__all__

# Then you can be explicit to control what ends up in the namespace,
# __all__ += ['do_primes']   # noqa
# or you can keep everything from the subpackage with the following instead
# __all__ += example_mod.__all__
