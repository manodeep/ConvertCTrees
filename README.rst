|Travis|

ConvertCTrees
#################

Converts Consistent Tree Output into LHaloTree format. *Written in C and probaly only works on Linux. Might depend on Linux kernel version*

Please **pay very close** attention to the Consistent-Trees version and **confirm** that the correct parsing is done in `here <https://github.com/manodeep/ConvertCTrees/blob/master/convert_trees_to_lhalo.c#L527>`_ (or `here <https://github.com/manodeep/ConvertCTrees/blob/master/convert_trees_to_lhalo.c#L473>`_, depending on the `Makefile` option `USE_STRINGPARSE`). Different versions of Consistent-Trees outputs different number of columns in the header file and the order assumed in the code might not be valid. See `Issue #4 <https://github.com/manodeep/ConvertCTrees/issues/4>`_



.. |Travis| image:: https://travis-ci.com/manodeep/ConvertCTrees.svg?branch=master
   :target: https://travis-ci.com/manodeep/ConvertCTrees
   :alt: Travis Compilation tests
