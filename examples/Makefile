# Build select examples for execution in a web environment

REMAT_INSTALL=../install/package
TARGETS=seismic_wave-web dissipative_wave-web potential_well-web contact-web buckling_arch-web bouncing_ball-web viscoplastic_truss-web triple_pendulum-web

all: $(TARGETS)

.PHONY: clean $(TARGETS)

$(TARGETS): %-web: %
# Make new sub-directory
	mkdir -p $@
# Copy any needed files from the base sub-directory
	cp -r $*/* $@/
# Rename the indicated example file as main.py
	mv $@/$*.py $@/main.py
# Copy the contents of the install/package directory to the sub-directory
	cp $(REMAT_INSTALL)/* $@
# Call pygbag to package project
	python3 -m pygbag --build $@

clean:
	rm -rf $(TARGETS)
