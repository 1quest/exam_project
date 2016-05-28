	#-------DEFINITIONS-----#
EXECUTABLE = ../opt_mini
INCLUDES = -lm

# ---- COMPILER AND FLAGS ---- #
CC := gcc
CFLAGS = -O3 -g -Wall
LDFLAGS = -lm
# ----------------------------- #

# ---- DIRECTORIES ---- #
UNOPDIR = u_opt
OPDIR = opt
BACK = ../
# --------------------- #

# ---- OBJECTS---- #
OBJS = funcs.o main.o ref_input.o
# ---------------- #

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(EXECUTABLE) $(INCLUDES)

$(OPDIR)/%.o: $(UNOPDIR)/%.c
	$(CC) $(CFLAGS) -c  $(LDFLAGS)

clean:
	rm -f *.o $(EXECUTABLE) $(EXECUTABLE).exe
