#pragma once

#include "PQPlanarityAttributes.h"

class PQPlanarity;

using namespace ogdf;

class PQPlanarityConsistency {
	PQPlanarity& pq;
	PQPlanarityDrawer draw;
	int checkCounter = 0;

public:
	static bool doWriteOut;

	explicit PQPlanarityConsistency(PQPlanarity& pq) : pq(pq), draw(&pq) {};

	bool consistencyCheck();

	void writeOut(std::string name = "", bool format = true, bool components = true);

	void checkComponentRegeneration();

	int getCheckCounter() const { return checkCounter; }
};
