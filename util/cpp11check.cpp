//*******************************************************************
// C++11 standard compliance check
//
// author: Carsten Gutwenger
//
// This test program checks if several C++11 features are supported
// by the C++ compiler. If the program compiles and executes without
// error message, the test is passed.
//
// The following C++11 features are checked:
//   (1)  Rvalue references (move semantics)
//   (2)  nullptr
//   (3)  Lambda expressions
//   (4)  Strongly typed enums
//   (5)  Right angle brackets
//   (6)  auto
//   (7)  Random number generators and distributions
//   (8)  Atomics
//   (9)  Time utilities
//   (10) Threads
//   (11) Range based for-loops
//   (12) Fixed-width integer types
//   (13) Initializer lists
//*******************************************************************

// C++ standard includes (must be available!)
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cmath>

#include <random>

#include <atomic>
#include <thread>
#include <chrono>

using std::cout;
using std::endl;
using std::flush;

bool passed = true;


void printCheck(bool ok)
{
	if(ok)
		cout << "ok." << endl;
	else {
		cout << "FAILED!" << endl;
		passed = false;
	}
}


// (1) Move semantics
template <class T> void myswap(T& a, T& b)
{
	T tmp(std::move(a));
	a = std::move(b);
	b = std::move(tmp);
}


// (3) Lambda expressions
void abssort(std::vector<float> &x)
{
	std::sort(x.begin(), x.end(),
		[](float a, float b) {
			return std::abs(a) < std::abs(b);
	});
}


// (6) auto

std::string hello("Hello auto!");

std::string &fAuto()
{
	return hello;
}


// (9) Time utilities
int fibonacci(int n)
{
	if (n < 3) return 1;
	return fibonacci(n-1) + fibonacci(n-2);
}


// (10) Threads
void computeSquares(int tid, std::vector<int> *squares)
{
	for(int i = tid*10; i < (tid+1)*10; ++i) {
		(*squares)[i] = i * i;
	}
}


// (13) Initializer lists
template<class T>
class MyContainer {
	T m_sum;

public:
	MyContainer(std::initializer_list<T> values) {
		m_sum = 0;
		for (auto x : values)
			m_sum += x;
	}

	T getSum() const { return m_sum; }
};


int main()
{
	//-----------------------------------------------------
	// (1) Check move semantics
	cout << "Check move semantics... " << flush;

	std::vector<int> a(100);
	std::vector<int> b(200);
	myswap(a,b);

	printCheck(a.size() == 200 && b.size() == 100);


	//-----------------------------------------------------
	// (2) Check nullptr
	cout << "Check nullptr and nullptr_t... " << flush;

	bool okNullptr = true;

	char* ch1 = nullptr;
	char* ch2 = nullptr;
	if(ch1 != ch2)
		okNullptr = false;
	if(ch1 != nullptr)
		okNullptr = false;
	if(ch1)
		okNullptr = false;

	std::nullptr_t nullVal = nullptr;
	char *ch3 = nullVal;
	if(ch3 != nullptr)
		okNullptr = false;

	printCheck(okNullptr);


	//-----------------------------------------------------
	// (3) Check lambda expressions
	cout << "Check lambda expressions... " << flush;

	bool okLambda = true;

	std::vector<float> vf(19);
	for(int i = 0; i < 19; ++i)
		vf[i] = float(i-9);

	abssort(vf);

	for(int i = 1; i < 19; ++i)
		if(std::abs(vf[i-1]) > std::abs(vf[i]))
			okLambda = false;

	printCheck(okLambda);


	//-----------------------------------------------------
	// (4) Check strongly typed enums
	cout << "Check strongly typed enums... " << flush;

	bool okEnums = true;

	enum class Colors { Green, Red, Blue, Berry = 100, Blueberry };
	enum class Days : short { Mo, Di, Mi, Do, Fr, Sa, So };

	if((int)Colors::Blueberry != 101)
		okEnums = false;

	enum E : unsigned long { E1 = 1, E2 = 2, Ebig = 0xFFFFFFF0U };
	if(Ebig != 0xFFFFFFF0U)
		okEnums = false;

	printCheck(okEnums);


	//-----------------------------------------------------
	// (5) Right angle brackets
	cout << "Check support for consecutive right angle brackets... " << flush;

	std::vector<std::vector<int>> table(4);
	table[2] = a;

	printCheck(table[2].size() == 200);


	//-----------------------------------------------------
	// (6) auto
	cout << "Check support for C++11 auto variable types... " << flush;

	bool okAuto = true;

	auto xAuto = 3.14;
	if(sizeof(xAuto) != sizeof(double))
		okAuto = false;

	auto        str1 = fAuto();
	auto       &str2 = fAuto();
	const auto &str3 = fAuto();

	if(str1 != hello || str2 != hello || str3 != hello)
		okAuto = false;

	printCheck(okAuto);


	//-----------------------------------------------------
	// (7) Random number generators and distributions
	cout << "Check random number generators and distributions... " << flush;

	bool okRand = true;

	std::minstd_rand rng;
	std::uniform_int_distribution<int> dist(1,6);

	std::vector<int> count(6);
	for(int i = 0; i < 100; ++i) {
		int r = dist(rng);
		if(r < 1 || r > 6)
			okRand = false;
		else
			count[r-1]++;
	}

#if 0
	for(int i = 0; i < 6; ++i)
		cout << " " << count[i];
#endif

	printCheck(okRand);


	//-----------------------------------------------------
	// (8) Atomics
	cout << "Check atomics... " << flush;

	bool okAtomic = true;

	std::atomic<int> atomicInt; // = 10;
	atomicInt = 10;
	atomicInt += 15;
	if(atomicInt != 25)
		okAtomic = false;

	int val = atomicInt.exchange(50);
	if(atomicInt != 50 || val != 25)
		okAtomic = false;

	printCheck(okAtomic);


	//-----------------------------------------------------
	// (9) Time utilities
	cout << "Check time utilities... " << flush;

	std::chrono::time_point<std::chrono::system_clock> start, end;

	start = std::chrono::system_clock::now();
	int result = fibonacci(36);
	end = std::chrono::system_clock::now();

	std::chrono::milliseconds::rep elapsed_milliseconds
		= std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
	cout << "[" << elapsed_milliseconds << " ms] ok" << endl;


	//-----------------------------------------------------
	// (10) Threads
	cout << "Check threads... " << flush;

	bool okThreads = true;

	const int numThreads = std::thread::hardware_concurrency();
	std::thread *t = new std::thread[numThreads];

	std::vector<int> squares(10*numThreads);
	for (int i = 0; i < numThreads; ++i) {
		t[i] = std::thread(computeSquares, i, &squares);
	}

	for (int i = 0; i < numThreads; ++i)
		t[i].join();

	delete [] t;

	for(std::vector<int>::size_type i = 0; i < squares.size(); ++i)
		if(squares[i] != i*i)
			okThreads = false;

	printCheck(okThreads);


	//-----------------------------------------------------
	// (11) Range based for-loops
	cout << "Check range based for-loops... " << flush;

	bool okRBFor = true;

	const int n = 10;
	std::vector<int> vec(n);
	for(int i = 0; i < n; ++i)
		vec[i] = 2*(i+1);

	int sum = 0;
	for(auto x : vec)
		sum += x;

	if(sum != n*(n+1))
		okRBFor = false;

	for(int &x : vec)
		x = 17;

	for(int x : vec)
		if(x != 17)
			okRBFor = false;

	printCheck(okRBFor);


	//-----------------------------------------------------
	// (12) Fixed-width integer types
	cout << "Check fixed-width integer types... " << flush;

	bool okFWInts = true;

	int8_t  xs8bit  = ( (int8_t ) 1 ) << 6;
	int16_t xs16bit = ( (int16_t) 1 ) << 14;
	int32_t xs32bit = ( (int32_t) 1 ) << 30;
	int64_t xs64bit = ( (int64_t) 1 ) << 62;

	okFWInts = okFWInts && (xs8bit > 0) && (xs16bit > 0) && (xs32bit > 0) && (xs64bit > 0);

	uint8_t  xus8bit  = ( (int8_t ) 1 ) << 7;
	uint16_t xus16bit = ( (int16_t) 1 ) << 15;
	uint32_t xus32bit = ( (int32_t) 1 ) << 31;
	uint64_t xus64bit = ( (int64_t) 1 ) << 63;

	okFWInts = okFWInts && (xus8bit > 0) && (xus16bit > 0) && (xus32bit > 0) && (xus64bit > 0);

	int_least8_t  xsl8bit  = ( (int8_t ) 1 ) << 6;
	int_least16_t xsl16bit = ( (int16_t) 1 ) << 14;
	int_least32_t xsl32bit = ( (int32_t) 1 ) << 30;
	int_least64_t xsl64bit = ( (int64_t) 1 ) << 62;

	okFWInts = okFWInts && (xsl8bit > 0) && (xsl16bit > 0) && (xsl32bit > 0) && (xsl64bit > 0);

	uint_least8_t  xusl8bit  = ( (int8_t ) 1 ) << 7;
	uint_least16_t xusl16bit = ( (int16_t) 1 ) << 15;
	uint_least32_t xusl32bit = ( (int32_t) 1 ) << 31;
	uint_least64_t xusl64bit = ( (int64_t) 1 ) << 63;

	okFWInts = okFWInts && (xusl8bit > 0) && (xusl16bit > 0) && (xusl32bit > 0) && (xusl64bit > 0);

	printCheck(okFWInts);


	//-----------------------------------------------------
	// (13) Initializer lists
	cout << "Check initializer lists for constructors... " << flush;

	bool okInitLists = true;

	MyContainer<int> myIntCont { 10, 11, 12, 13, 14, 15, 16, 17 };
	if (myIntCont.getSum() != 10 + 11 + 12 + 13 + 14 + 15 + 16 + 17)
		okInitLists = false;

	printCheck(okInitLists);


	//-----------------------------------------------------
	cout << "\nCompliance check ";
	if(passed)
		cout << "passed.";
	else
		cout << "FAILED!";
	cout << endl;


	return 0;
}
