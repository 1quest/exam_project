#-------DEFINITIONS-----#
EXECUTABLE = ../unopt_mini
INCLUDES = -lm

# ---- COMPILER AND FLAGS ---- #
CC := gcc
CFLAGS = -g -Wall
LDFLAGS = -lm
# ----------------------------- #

# ---- DIRECTORIES ---- #
UNOPDIR = u_opt
OPDIR = opt
BACK = ../
# --------------------- #

# ---- OBJECTS---- #
OBJS = funcs.o main.o
# ---------------- #

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(EXECUTABLE) $(INCLUDES)

$(UNOPDIR)/%.o: $(UNOPDIR)/%.c
	$(CC) $(CFLAGS) -c  $(LDFLAGS)

clean:
	rm -f *.o $(EXECUTABLE)
