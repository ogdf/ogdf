all: all-basic all-layout all-special

clean: clean-basic clean-layout clean-special

all-%:
	$(MAKE) -C $(patsubst all-%,%,$@) all

clean-%:
	$(MAKE) -C $(patsubst clean-%,%,$@) clean
