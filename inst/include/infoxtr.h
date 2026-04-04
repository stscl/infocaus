// Copyright (C) 2026 Wenbo Lyu
//
// This file is part of infoxtr.
//
// infoxtr is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// infoxtr is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with infoxtr. If not, see <https://www.gnu.org/licenses/>.

#ifndef INFOXTR_INFOXTR_H
#define INFOXTR_INFOXTR_H

// ============================================================
// Dependency Guard: Encourage best practices (non-blocking)
// ============================================================
#if defined(Rcpp_hpp) && !defined(COMPILING_INFOXTR)
    #warning "It is recommended to include <infoxtr.h> alone, as it already includes <Rcpp.h> and <RcppThread.h>."
#endif

// ============================================================
// Core Dependencies (Auto-included for users)
// ============================================================

#include <Rcpp.h>
#include <RcppThread.h>

// ============================================================
// Module Headers (Organized by functionality)
// ============================================================

#include "infoxtr/lagg.hpp"
#include "infoxtr/embed.hpp"
#include "infoxtr/combn.hpp"
#include "infoxtr/numericutils.hpp"
#include "infoxtr/distance.hpp"
#include "infoxtr/neighbor.hpp"
#include "infoxtr/discretize.hpp"
#include "infoxtr/infotheo.hpp"
#include "infoxtr/ksginfo.hpp"
#include "infoxtr/transferentropy.hpp"
#include "infoxtr/surd.hpp"

#endif // INFOXTR_INFOXTR_H