
all-local:
	cd Cuba-$(CUBA_VERSION) && $(MAKE) $(AM_MAKEFLAGS) lib
check-local:
#	cd Cuba-$(CUBA_VERSION) && $(MAKE) $(AM_MAKEFLAGS) test
clean-local:
	cd Cuba-$(CUBA_VERSION) && $(MAKE) $(AM_MAKEFLAGS) clean

# Assuming the package knows how to install itself
install-data-local:
#	cd Cuba-$(CUBA_VERSION) && $(MAKE) $(AM_MAKEFLAGS) install-data
install-exec-local:
#	cd Cuba-$(CUBA_VERSION) && $(MAKE) $(AM_MAKEFLAGS) install
uninstall-local:
#	cd Cuba-$(CUBA_VERSION) && $(MAKE) $(AM_MAKEFLAGS) uninstall

distclean-local:
	rm -rf Cuba-4.2/config.status Cuba-4.2/config.log Cuba-4.2/libcuba.a Cuba-4.2/makefile Cuba-4.2/config.h

# Distribute files from here.
EXTRA_DIST = Cuba-$(CUBA_VERSION)/src Cuba-$(CUBA_VERSION)/tools Cuba-$(CUBA_VERSION)/configure.ac Cuba-$(CUBA_VERSION)/configure Cuba-$(CUBA_VERSION)/makefile.in Cuba-$(CUBA_VERSION)/config.h.in 

pkginclude_HEADERS = Cuba-$(CUBA_VERSION)/cuba.h

