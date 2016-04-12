#OPT += -DUSE_FGETS
#OPT += -USE_MPI

include common.mk

target := convert_trees_to_lhalo
SRC1   := convert_trees_to_lhalo.c utils.c progressbar.c 
OBJS1  := $(SRC1:.c=.o)
INCL   := output_datatype.h utils.h progressbar.h 

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


