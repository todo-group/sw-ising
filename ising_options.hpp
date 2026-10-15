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

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>

struct options {
  unsigned int seed, q, length_x, length_y;
  double temperature;
  double coupling[2] = {1.0, 1.0};
  unsigned int sweeps, therm;
  bool valid, print_graph;

  options(unsigned int argc, char *argv[], bool print = true) :
    // default parameters
    seed(29833), length_x(4), length_y(4), temperature(2.27), sweeps(1 << 16), therm(sweeps >> 3),
    valid(true), print_graph(false) {
    for (unsigned int i = 1; i < argc; ++i) {
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 's' :
          if (++i == argc) { usage(print); return; }
          seed = std::atoi(argv[i]); break;
        case 'l' : {
          unsigned int* length;
          if (std::strcmp(argv[i], "-lx") == 0) length = &length_x;
          else if (std::strcmp(argv[i], "-ly") == 0) length = &length_y;
          else { usage(print); return; }
          if (++i == argc) { usage(print); return; }
          *length = std::atoi(argv[i]); break;
        }
        case 't' :
          if (++i == argc) { usage(print); return; }
          temperature = std::atof(argv[i]); break;
        case 'j' : {
          int type;
          if (std::strcmp(argv[i], "-jx") == 0) type = 0;
          else if (std::strcmp(argv[i], "-jy") == 0) type = 1;
          else { usage(print); return; }
          if (++i == argc) { usage(print); return; }
          char* end;
          coupling[type] = std::strtod(argv[i], &end);
          if (end == argv[i] || *end != '\0') { usage(print); return; }
          break;
        }
        case 'm' :
          if (++i == argc) { usage(print); return; }
          sweeps = std::atoi(argv[i]);
          therm = sweeps >> 3; break;
        case 'g' :
          print_graph = true; break;
        case 'h' :
          usage(print, std::cout); return;
        default :
          usage(print); return;
        }
        break;
      default :
        usage(print); return;
      }
    }
    if (length_x == 0 || length_y == 0 || temperature <= 0. || sweeps == 0 ||
        !std::isfinite(coupling[0]) || !std::isfinite(coupling[1]) ||
        coupling[0] < 0. || coupling[1] < 0.) {
      std::cerr << "invalid parameter(s)\n"; usage(print); return;
    }
    if (print) {
      std::cout << "Seed of RNG            = " << seed << std::endl
                << "System Size (Lx)       = " << length_x << std::endl
                << "System Size (Ly)       = " << length_y << std::endl
                << "Temperature            = " << temperature << std::endl
                << "Coupling Jx            = " << coupling[0] << std::endl
                << "Coupling Jy            = " << coupling[1] << std::endl
                << "MCS for Thermalization = " << therm << std::endl
                << "MCS for Measurement    = " << sweeps << std::endl;
    }
  }
  void usage(bool print, std::ostream& os = std::cerr) {
    if (print)
      os << "[command line options]\n"
         << "  -s int    Seed of RNG\n"
         << "  -lx int   System Size Lx (default 4)\n"
         << "  -ly int   System Size Ly (default 4)\n"
         << "  -t double Temperature\n"
         << "  -jx double Coupling for bond type 0 (x, nonnegative; default 1.0)\n"
         << "  -jy double Coupling for bond type 1 (y, nonnegative; default 1.0)\n"
         << "  -m int    MCS for Measurement\n"
         << "  -g        Print lattice information and exit\n"
         << "  -h        this help\n";
    valid = false;
  }
};
