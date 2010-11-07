all:
	$(MAKE) -C src
clean:
	cd src && $(MAKE) clean
