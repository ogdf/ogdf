/** \file
 *
 * \author Marcel Radermacher
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

#include <ogdf/basic/basic.h>

#include <cstddef>
#include <limits>
#include <vector>

namespace ogdf {
namespace internal {
namespace gcm {
namespace datastructure {

class TimestampFlags {
private:
	std::vector<unsigned int> flags;
	unsigned int current_round = 1;

public:
	TimestampFlags() { /*nothing to do*/
	}

	TimestampFlags(size_t size) : flags(size, 0) { /*nothing to do*/
	}

	size_t size() const { return flags.size(); }

	void clear() {
		++current_round;
		if (current_round == std::numeric_limits<unsigned int>::max()) {
			current_round = 1;
			std::fill(flags.begin(), flags.end(), 0);
		}
	}

	inline bool is_set(size_t id) const {
		OGDF_ASSERT(id < flags.size());
		return flags[id] == current_round;
	}

	inline void set(size_t id) {
		OGDF_ASSERT(id < flags.size());
		flags[id] = current_round;
	}
};

}
}
}
}
