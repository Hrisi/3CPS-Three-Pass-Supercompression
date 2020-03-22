all:
	@+make --no-print-directory -C build-debug/
install:
	@+make install --no-print-directory -C build-debug/
package:
	@+make package --no-print-directory -C build-debug/
test:
	@+make test --no-print-directory -C build-debug/
clean:
	@+make clean --no-print-directory -C build-debug/
distclean:
	@rm -Rf build-debug/
