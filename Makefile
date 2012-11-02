CC=mpicc
CFLAGS=-c
LDFLAGS=-lm
EXEC=nbody
SRCS=nbody.c
OBJS=$(SRCS:.c=.o)
$(EXEC):$(OBJS)
	$(CC) $^ -o $@ $(LDFLAGS)
.c.o:
	$(CC) $(CFLAGS) $< -o $@
