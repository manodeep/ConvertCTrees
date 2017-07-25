# ConvertCTrees

# Description 

Converts Consistent Tree Output into LHaloTree format. *Written in C and probaly only works on Linux. Might depend on Linux kernel version*

Please **pay very close** attention to the Consistent-Trees version and **confirm** that the correct parsing is done in [here](https://github.com/manodeep/ConvertCTrees/blob/master/convert_trees_to_lhalo.c#L527) (or [here](https://github.com/manodeep/ConvertCTrees/blob/master/convert_trees_to_lhalo.c#L473), depending on the `Makefile` option `USE_STRINGPARSE`). Different versions of Consistent-Trees outputs different number of columns in the header file and the order assumed in the code might not be valid. See [here](https://github.com/manodeep/ConvertCTrees/issues/4)



