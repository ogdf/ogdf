include(CheckCXXSourceCompiles)

check_cxx_source_compiles("
#include <malloc.h>
int main() {
	struct mallinfo2 mi = mallinfo2();
	return 0;
}" has_mallinfo2)
