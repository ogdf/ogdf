/** \file
 * \brief TODO Document
 *
 * \author Simon D. Fink <ogdf@niko.fink.bayern>
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.md in the OGDF root directory for details.
 *
 * \par
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * Version 2 or 3 as published by the Free Software Foundation;
 * see the file LICENSE.txt included in the packaging of this file
 * for details.
 *
 * \par
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, see
 * http://www.gnu.org/copyleft/gpl.html
 */
#pragma once

#include <json.hpp>

using nlohmann::json;

#ifndef LIKWID_PERFMON

void likwid_prepare(json& results) { }

void likwid_finalize(json& results) { }

#else

#	include <iostream>
#	include <sstream>
#	include <string>
#	include <vector>

#	include <likwid.h>

#	define MAX_NUM_EVENTS 10

struct RegionData {
	std::vector<double> events;
	double time = 0;
	int count = 0;

	void get(const std::string& region) {
		int nevents = MAX_NUM_EVENTS;
		double events_a[MAX_NUM_EVENTS] = {0.0};
		likwid_markerGetRegion(region.c_str(), &nevents, events_a, &time, &count);
		events.assign(events_a, events_a + nevents);
		// printf("Region '%s' was called %d times, measured %d events, total measurement time is %f\n", region.c_str(), count, nevents, time);
	}

	void diff_update(const RegionData& other) {
		//OGDF_ASSERT(events.size() == other.events.size());
		time = other.time - time;
		count = other.count - count;
		for (int i = 0; i < events.size(); i++) {
			events[i] = other.events[i] - events[i];
		}
	}

	json toJson() const {
		json ret;
		ret["nevents"] = events.size();
		ret["time"] = time;
		ret["count"] = count;
		ret["data"] = events;
		return ret;
	}

	json toJson(std::vector<std::string>& event_names) const {
		//OGDF_ASSERT(event_names.size() == events.size());
		json ret;
		ret["nevents"] = events.size();
		ret["time"] = time;
		ret["count"] = count;
		json& data = ret["data"] = json::object();
		for (int i = 0; i < events.size(); i++) {
			data[event_names[i]] = events[i];
		}
		return ret;
	}
};

std::vector<std::string> regions {"makeReduced", "makeReduced-step", "solveReduced", "embed",
		"embed-step",

		"convertSmall", "makeWheel", "contractWheel", "contract", "contract-bicon",
		"contract-encapsulate", "encapsulate", "checkPCTree", "propagatePQ", "propagatePQ-makePCv",
		"simplify", "simplify-bondmap", "simplify-trans", "batchSPQR", "batchSPQR-makeTree",
		"batchSPQR-makeSPQR", "batchSPQR-simplify", "batchSPQR-propagatePQ",

		"undo-convertSmall", "undo-contract", "undo-contract-collect-bicon",
		"undo-contract-collect-stars", "undo-encapsulate", "undo-propagatePQ", "undo-simplify",

		"ogdf-array-construct", "ogdf-array-fill", "ogdf-array-deconstruct"};
std::unordered_map<std::string, RegionData> likwidRegions;
std::vector<std::string> likwidEvents, likwidCounters, likwidMetrics;

void likwid_prepare(json& results) {
	results["likwid_perfmon"] = true;
	likwid_markerInit();
	likwid_markerThreadInit();
	timer_init();
	for (std::string& region : regions) {
		likwid_markerRegisterRegion(region.c_str());
		likwidRegions[region].get(region);
	}
	int nevents = perfmon_getNumberOfEvents(perfmon_getIdOfActiveGroup());
	likwidEvents.reserve(nevents);
	likwidCounters.reserve(nevents);
	likwidMetrics.reserve(nevents);
	for (int i = 0; i < nevents; i++) {
		likwidEvents.emplace_back(perfmon_getEventName(perfmon_getIdOfActiveGroup(), i));
		likwidCounters.emplace_back(perfmon_getCounterName(perfmon_getIdOfActiveGroup(), i));
	}
	for (int i = 0; i < perfmon_getNumberOfMetrics(perfmon_getIdOfActiveGroup()); i++) {
		likwidMetrics.emplace_back(perfmon_getMetricName(perfmon_getIdOfActiveGroup(), i));
	}
	results["likwid_meta"] = {
			{"cpu_freq", timer_getCpuClock()},
			{"group_id", perfmon_getIdOfActiveGroup()},
			{"ngroups", perfmon_getNumberOfGroups()},
			{"nthreads", perfmon_getNumberOfThreads()},
			{"nregions", perfmon_getNumberOfRegions()},
			{"nmetrics", perfmon_getNumberOfMetrics(perfmon_getIdOfActiveGroup())},
			{"group_name", perfmon_getGroupName(perfmon_getIdOfActiveGroup())},
			{"group_info_s", perfmon_getGroupInfoShort(perfmon_getIdOfActiveGroup())},
			{"group_info_l", perfmon_getGroupInfoLong(perfmon_getIdOfActiveGroup())},
			{"nevents", nevents},
			{"event_names", likwidEvents},
			{"metric_names", likwidMetrics},
			{"counter_names", likwidCounters},
	};
}

void likwid_finalize(json& results) {
	json& regionsJson = results["likwid"] = json::object();
	for (std::string& region : regions) {
		likwidRegions[region].get(region);
		regionsJson[region] = likwidRegions[region].toJson();
	}
	timer_finalize();
	likwid_markerClose();
}

#endif
