# Build select examples for execution in a web environment

REMAT_INSTALL=../install/package
TARGETS=double_pendulum-web # contact-web

all: $(TARGETS)

.PHONY: clean $(TARGETS)

$(TARGETS): %-web: %.py
# Make new directory
	mkdir -p $@
# Copy the indicated example file into the sub-directory and rename as main.py
	cp $< $@/main.py
# Copy the contents of the install/package directory to the sub-directory
	cp $(REMAT_INSTALL)/* $@
# Call pygbag to package project
	python3 -m pygbag --build $@

clean:
	rm -rf $(TARGETS)
