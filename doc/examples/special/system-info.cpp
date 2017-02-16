#include <ogdf/basic/basic.h>
#include <ogdf/basic/System.h>

using namespace ogdf;

const char *yn(bool b)
{
	return (b ? "yes" : "no");
}

int main()
{
	cout << "---------------------------------------" << endl;
	cout << "      System-specific information      " << endl;
	cout << "---------------------------------------" << endl;
	cout << endl;

	cout << "Cache / processors:" << endl;
	cout << "-------------------" << endl;
	cout << "Processors: " << System::numberOfProcessors() << endl;
	cout << "L2-Cache:   " << System::cacheSizeKBytes() << " KBytes" << endl;
	cout << "Cache-Line: " << System::cacheLineBytes() << " Bytes" << endl;
	cout << endl;

	cout << "Supported technologies:" << endl;
	cout << "-----------------------" << endl;
	cout << "MMX:    " << yn(System::cpuSupports(CPUFeature::MMX))    << endl;
	cout << "SSE:    " << yn(System::cpuSupports(CPUFeature::SSE))    << endl;
	cout << "SSE2:   " << yn(System::cpuSupports(CPUFeature::SSE2))   << endl;
	cout << "SSE3:   " << yn(System::cpuSupports(CPUFeature::SSE3))   << endl;
	cout << "SSSE3:  " << yn(System::cpuSupports(CPUFeature::SSSE3))  << endl;
	cout << "SSE4.1: " << yn(System::cpuSupports(CPUFeature::SSE4_1)) << endl;
	cout << "SSE4.2: " << yn(System::cpuSupports(CPUFeature::SSE4_2)) << endl;
	cout << "VMX:    " << yn(System::cpuSupports(CPUFeature::VMX))    << endl;
	cout << "SMX:    " << yn(System::cpuSupports(CPUFeature::SMX))    << endl;
	cout << "EST:    " << yn(System::cpuSupports(CPUFeature::EST))    << endl;
	cout << endl;

	cout << "Memory management:" << endl;
	cout << "------------------" << endl;
	cout << "Total physical memory: " << System::physicalMemory() / 1024 / 1024 << " MBytes" << endl;
	cout << "  available:           " << System::availablePhysicalMemory() / 1024 / 1024 << " MBytes" << endl;
	cout << "  used by process:     " << System::memoryUsedByProcess() / 1024 << " KBytes" << endl;
#if defined(OGDF_SYSTEM_WINDOWS) || defined(__CYGWIN__)
	cout << "  peak amount:         " << System::peakMemoryUsedByProcess() / 1024 << " KBytes" << endl;
#endif
	cout << endl;
	cout << "allocated by malloc:   " << System::memoryAllocatedByMalloc() / 1024 << " KBytes" << endl;
	cout << "  in freelist:         " << System::memoryInFreelistOfMalloc() / 1024 << " KBytes" << endl;
	cout << endl;
	cout << "allocated by OGDF:     " << System::memoryAllocatedByMemoryManager() / 1024 << " KBytes" << endl;
	cout << "  in global freelist:  " << System::memoryInGlobalFreeListOfMemoryManager() / 1024 << " KBytes" << endl;
	cout << "  in thread freelist:  " << System::memoryInThreadFreeListOfMemoryManager() / 1024 << " KBytes" << endl;

	return 0;
}
