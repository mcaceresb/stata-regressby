# ---------------------------------------------------------------------
# OS parsing

ifeq ($(OS),Windows_NT)
	OSFLAGS = -shared
	GCC = x86_64-w64-mingw32-gcc.exe
	OUT = regressby_windows.plugin
	GSL = -L/lib -l:libgsl.dll.a -l:libgslcblas.dll.a
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		OSFLAGS = -shared -fPIC -DSYSTEM=OPUNIX
		OUT = regressby_unix.plugin
	endif
	ifeq ($(UNAME_S),Darwin)
		OSFLAGS = -bundle -DSYSTEM=APPLEMAC
		OUT = regressby_macosx.plugin
	endif
	GSL = -L/usr/local/lib -l:libgsl.a -l:libgslcblas.a
	GCC = gcc
endif

CFLAGS = -Wall -O3 $(OSFLAGS)

# ---------------------------------------------------------------------
# Main

all: clean regressby

regressby: plugin/stplugin.c plugin/regressby.c
	$(GCC) $(CFLAGS) -o $(OUT) plugin/stplugin.c plugin/regressby.c $(GSL)

.PHONY: clean
clean:
	rm -f $(OUT)
