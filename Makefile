# ---- DIRECTORIES ---- #
UNOPDIR = u_opt
OPDIR = opt
BACK = ../
# --------------------- #

all: $(UNOPDIR)

$(UNOPDIR):
	$(MAKE) -C $(UNOPDIR)

.PHONY: $(UNOPDIR)

clean:
	$(MAKE) clean -C $(UNOPDIR)
