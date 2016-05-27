# ---- DIRECTORIES ---- #
UNOPDIR = u_opt
OPDIR = opt
BACK = ../
# --------------------- #

all: $(UNOPDIR)

$(UNOPDIR):
	$(MAKE) -C $(UNOPDIR)
	$(MAKE) -C $(OPDIR)

.PHONY: $(UNOPDIR)
.PHONY: $(OPDIR)

clean:
	$(MAKE) clean -C $(UNOPDIR)
	$(MAKE) clean -C $(OPDIR)
