/*****************************************************************************
*
* Cluster-MC: Cluster Algorithm Monte Carlo Methods
*
* Copyright (C) 1997-2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     https://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*
*****************************************************************************/

// Swendsen-Wang Algorithm for 2D square-lattice Ising model with anisotropic interaction

#ifndef ALPS_INDEP_SOURCE
# define ALPS_INDEP_SOURCE
#endif

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <standards/accumulator.hpp>
#include <standards/power.hpp>
#include <standards/timer.hpp>
#include <lattis/graph.hpp>
#include "union_find.hpp"
#include "ising_options.hpp"

using standards::power2;
using standards::power4;

int main(int argc, char* argv[]) {
  std::cout << "Swendsen-Wang Algorithm for 2D square-lattice Ising model with anisotropic interaction\n";
  options p(argc, argv);
  if (!p.valid) std::exit(127);
  double prob[2] = {1 - std::exp(-2 * p.coupling[0] / p.temperature),
                    1 - std::exp(-2 * p.coupling[1] / p.temperature)};

  // rectangular lattice: x bonds have type 0, y bonds have type 1
  lattis::unitcell cell(2);
  cell.add_site(lattis::coordinate_t::Zero(2), 0);
  for (int direction = 0; direction < 2; ++direction) {
    lattis::offset_t offset = lattis::offset_t::Zero(2);
    offset(direction) = 1;
    cell.add_bond(0, 0, offset, direction);
  }
  lattis::span_t span = lattis::span_t::Zero(2, 2);
  span(0, 0) = p.length_x;
  span(1, 1) = p.length_y;
  lattis::graph lattice(lattis::basis::simple(2), cell, span);
  if (p.print_graph) {
    lattice.print();
    for (std::size_t b = 0; b < lattice.num_bonds(); ++b) {
      std::cout << "bond: " << b << " type: " << lattice.bond_type(b)
                << " source: " << lattice.source(b)
                << " target: " << lattice.target(b) << '\n';
    }
    return 0;
  }

  // random number generators
  std::mt19937 eng(p.seed);
  std::uniform_real_distribution<> r_uniform01;

  // spin configuration
  std::vector<int> spins(lattice.num_sites(), 1);

  // cluster information
  typedef cluster::union_find::node fragment_t;
  std::vector<fragment_t> fragments(lattice.num_sites());
  std::vector<int> flip(lattice.num_sites());

  // observables
  standards::accumulator num_clusters("Number of Clusters"), energy("Energy Density"),
    magnetization_unimp("Magnetization (unimproved)"),
    magnetization2_unimp("Magnetization^2 (unimproved)"),
    magnetization4_unimp("Magnetization^4 (unimproved)"),
    magnetization2("Magnetization^2"), magnetization4("Magnetization^4");

  standards::timer tm;
  for (unsigned int mcs = 0; mcs < p.therm + p.sweeps; ++mcs) {
    // initialize cluster information
    std::fill(fragments.begin(), fragments.end(), fragment_t());

    // cluster generation
    for (int b = 0; b < lattice.num_bonds(); ++b) {
      if (spins[lattice.source(b)] == spins[lattice.target(b)] &&
          r_uniform01(eng) < prob[lattice.bond_type(b)])
        unify(fragments, lattice.source(b), lattice.target(b));
    }

    // assign cluster id & accumulate cluster properties
    int nc = 0;
    double mag2 = 0, mag4 = 0;
    for (auto& f : fragments) {
      if (f.is_root()) {
        f.set_id(nc++);
        double w = f.weight();
        mag2 += power2(w);
        mag4 += power4(w);
      }
    }
    for (auto& f : fragments) f.set_id(cluster_id(fragments, f));

    // flip spins
    for (int c = 0; c < nc; ++c) flip[c] = (r_uniform01(eng) < 0.5);
    for (int s = 0; s < lattice.num_sites(); ++s)
      if (flip[fragments[s].id()]) spins[s] ^= 1;

    double ene = 0;
    for (int b = 0; b < lattice.num_bonds(); ++b) {
      ene -= p.coupling[lattice.bond_type(b)] *
        (spins[lattice.source(b)] == spins[lattice.target(b)] ? 1.0 : -1.0);
    }
    double mu = 0;
    for (int s = 0; s < lattice.num_sites(); ++s) mu += 2 * spins[s] - 1;

    if (mcs >= p.therm) {
      num_clusters << (double)nc;
      energy << ene / lattice.num_sites();
      magnetization_unimp << mu;
      magnetization2_unimp << power2(mu);
      magnetization4_unimp << power4(mu);
      magnetization2 << mag2;
      magnetization4 << (3 * power2(mag2) - 2 * mag4);
    }
  }

  double elapsed = tm.elapsed();
  std::clog << "Elapsed time = " << elapsed << " sec\n"
            << "Speed = " << (p.therm + p.sweeps) / elapsed << " MCS/sec\n";
  std::cout << num_clusters << std::endl
            << energy << std::endl
            << magnetization_unimp << std::endl
            << magnetization2_unimp << std::endl
            << magnetization4_unimp << std::endl
            << magnetization2 << std::endl
            << magnetization4 << std::endl
            << "Binder Ratio of Magnetization = "
            << power2(magnetization2.mean()) / magnetization4.mean() << std::endl;
}
