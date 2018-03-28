#include <ogdf/basic/basic.h>

int main() {
#ifdef OGDF_DEBUG
	bool debugMode = true;
#else
	bool debugMode = false;
#endif
	std::cout << "This user program is compiled in " << (debugMode ? "Debug" : "Release") << " mode." << std::endl;
	std::cout << "The OGDF is compiled in " << (ogdf::debugMode ? "Debug" : "Release") << " mode." << std::endl;
	if (debugMode != ogdf::debugMode) {
		std::cout << "Check your configuration!" << std::endl;
		return 1;
	}
	std::cout << "Everything is fine!" << std::endl;
	return 0;
}
