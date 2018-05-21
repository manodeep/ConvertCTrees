## These are code optimization options
#OPT += -DUSE_FGETS
#OPT += -USE_MPI
OPT += -DUSE_STRINGPARSE


# This is to correct the snapshot index output by
# Consistent Trees into the actual simulation snapshot
# number. See https://bitbucket.org/pbehroozi/consistent-trees/commits/41866601a537ce5f3d159d0364ef7d7082da43a5
# For instance, if the first 5 snapshots in the simulation
# did not contain sufficient halos, then the Consistent Trees `snap_idx == 0`
# will correspond to the simulation `snapnum == 5`
OPT += -DSNAP_OFFSET=0

include common.mk

target := convert_trees_to_lhalo
SRC1   := main.c convert_trees_to_lhalo.c utils.c progressbar.c
INCL   := output_datatype.h utils.h progressbar.h stringparse.h check_syscalls.h convert_trees_to_lhalo.h
ifeq (USE_STRINGPARSE,$(findstring USE_STRINGPARSE,$(OPT)))
  SRC1 += stringparse.c check_syscalls.c
  INCL += stringparse.h check_syscalls.h strtonum.c
endif
OBJS1  := $(SRC1:.c=.o)

SRC2  := tests/compare_two_lhalotree.c utils.c 
OBJS2 := $(SRC2:.c=.o)

all: $(target) $(SRC1) $(INCL) test Makefile 

$(target): $(OBJS1) $(INCL) Makefile common.mk
	$(CC) $(OBJS1) $(CLINK) $(CFLAGS) -o $@

%.o: %.c $(INCL) Makefile common.mk
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

test: compare_two_lhalotree $(INCL) common.mk Makefile

compare_two_lhalotree: $(OBJS2) $(INCL) common.mk Makefile
	$(CC) $(OBJS2) $(CLINK) $(CFLAGS) -o $@

.PHONY: clean clena celan 

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2)

clena: clean

celan: celan



