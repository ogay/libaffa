lib_LTLIBRARIES = libaa.la
libaa_la_SOURCES = aa_rounding.cpp aa_interval.cpp aa_aaftrigo.cpp aa_aafapprox.cpp aa_aafarithm.cpp aa_aafcommon.cpp
libaa_la_LDFLAGS =      \
	-release $(LT_RELEASE)  \
        -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)

include_HEADERS = aa.h aa_aaf.h aa_interval.h
noinst_HEADERS = aa_rounding.h
