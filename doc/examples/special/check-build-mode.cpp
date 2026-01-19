#include <ogdf/basic/basic.h>
#include <iostream>

int main() {
#ifdef OGDF_DEBUG
	bool debugMode = true;
#else
	bool debugMode = false;
#endif
	std::cout << "This user program uses includes/headers for OGDF " << (debugMode ? "Debug" : "Release") << " mode." << std::endl;
	std::cout << "The linked OGDF library was compiled in " << (ogdf::debugMode ? "Debug" : "Release") << " mode." << std::endl;
	if (debugMode != ogdf::debugMode) {
		std::cout << "Check your configuration!" << std::endl;
		return 1;
	}
	std::cout << "Everything is fine!" << std::endl;
	return 0;
}
