# from __future__ import absolute_import
from __future__ import division # confidence high

from .version import *

try:
    from reftools.svninfo import (__svn_version__, __full_svn_info__,
                                  __setup_datetime__)
except ImportError:
    pass
