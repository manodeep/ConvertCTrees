#OPT += -DUSE_FGETS
#OPT += -USE_MPI

include common.mk

target=convert_trees_to_lhalo
SRC1=convert_trees_to_lhalo.c utils.c progressbar.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = output_datatype.h utils.h progressbar.h 

all: $(target) $(SRC1) $(INCL) Makefile 

$(target): $(OBJS1) $(INCL) Makefile common.mk
	$(CC) $(OBJS1) $(CLINK) $(CFLAGS) -o $@

%.o: %.c $(INCL) Makefile common.mk
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean clena celan 

clean:
	$(RM) $(target) $(OBJS1) 

clena: clean

celan: celan



