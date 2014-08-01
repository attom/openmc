all:
	mkdir -p build
	cmake -H. -Bbuild
	make -s -C build
clean:
	make -s -C build clean
distclean:
	rm -fr build
test:
	make -s -C build test
install:
	make -s -C build install
till:
	gfortran clustering.F90 -o clustering -fbounds-check

.PHONY: all clean distclean test install
