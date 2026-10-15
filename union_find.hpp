/*****************************************************************************
*
* Cluster-MC: Cluster Algorithm Monte Carlo Methods
*
* Copyright (C) 1997-2019 by Synge Todo <wistaria@comp-phys.org>
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

// Thread-safe weighted union-find algorithm
// Reference:
//   D. Knuth,
//   `The Art of Computer Programming, Vol. 1, Fundamental Algorithms'
//   3rd edition (Addison Wesley, Reading, 1997) Sec 2.3.3.

#pragma once

#include <algorithm> // for std::swap
#include <iostream>
#include <limits>
#include <vector>

#ifdef _OPENMP
# include <omp.h>
# include <atomic>
#endif

namespace cluster {
namespace union_find {

class node {
public:
  node() : parent_(-1) {} // root node with weight = 1
#ifdef _OPENMP
  node(const node& rhs) : id_(rhs.id_) { parent_.store(rhs.parent_.load()); }
  node& operator=(const node& rhs) {
    parent_.store(rhs.parent_.load());
    id_ = rhs.id_;
    return *this;
  }
#else
  // default copy constructor and assignment operator is fine
#endif
  bool is_root() const { return parent_ <= 0; }
  void set_parent(int parent) { parent_ = parent + 1; }
  int parent() const { return parent_ - 1; }
  void set_weight(int w) { parent_ = -w; }
  int weight() const { return -parent_; }
  void set_id(int id) { id_ = id; }
  int id() const { return id_; }
#ifdef _OPENMP
  int lock_root() {
    int p = parent_;
    return (p < 0 && parent_.compare_exchange_weak(p, 0)) ? (-p) : 0;
    //// example of THREAD-UNSAFE implementation
    //   if (parent_ < 0) {
    //     int p = parent_;
    //     parent_ = 0;
    //     return -p;
    //   } else {
    //     return 0;
    //   }
  }
  // unlock can be done by set_parent or set_weight
#endif
private:
#ifdef _OPENMP
  std::atomic<int> parent_; // negative for root fragment, zero for locked root
#else
  int parent_;
#endif
  int id_;
};

class node_noweight {
public:
  node_noweight() : parent_(id_mask) {} // root note with id = 0
#ifdef _OPENMP
  node_noweight(const node_noweight& rhs) { parent_.store(rhs.parent_.load()); }
  node_noweight& operator=(const node_noweight& rhs) {
    parent_.store(rhs.parent_.load());
    return *this;
  }
#else
  // default copy constructor and assignment operator is fine
#endif
  int is_root() const { return (int) parent_ <= 0; }
  void set_parent(int parent) { parent_ = parent + 1; }
  int parent() const { return parent_ - 1; }
  void set_weight(int) { set_id(0); } // dummy routine for unlock
  int weight() const { return 0; } // dummy
  void set_id(int id) { parent_ = id ^ id_mask; }
  int id() const { return parent_ ^ id_mask; }
#ifdef _OPENMP
  int lock_root() {
    int p = parent_;
    return (p < 0 && parent_.compare_exchange_weak(p, 0)) ? (-p) : 0;
  }
  // unlock can be done by set_parent, set_weight or set_id
#endif
private:
  static const int id_mask = -1;
#ifdef _OPENMP
  std::atomic<int> parent_; // negative for root fragment, zero for locked root
#else
  int parent_;
#endif
};

// thread-unsafe
template<class T>
inline int add(std::vector<T>& v) {
  v.push_back(T());
  return v.size() - 1; // return index of new node
}

template<class T>
inline int root_index(std::vector<T> const& v, int g) {
#ifdef _OPENMP
  T c = v[g];
  while (!c.is_root()) {
    g = c.parent();
    c = v[g];
  }
#else
  while (!v[g].is_root()) {
    g = v[g].parent();
  }
#endif
  return g;
}

// root_index with path-halving
// Note: this is not thread-safe, but is really safe as long as called from unify_*
template<class T>
inline int root_index_ph(std::vector<T>& v, int g) {
  if (v[g].is_root()) return g;
  while (true) {
    int p = v[g].parent();
    if (v[p].is_root()) return p;
    v[g].set_parent(v[p].parent());
    g = p;
  }
}

template<class T>
inline T const& root(std::vector<T> const& v, int g) { return v[root_index(v, g)]; }

template<class T>
inline T const& root(std::vector<T> const& v, T const& n) {
  return n.is_root() ? n : root(v, n.parent());
}

template<class T>
inline int cluster_id(std::vector<T> const& v, int g) { return root(v, g).id(); }

template<class T>
inline int cluster_id(std::vector<T> const& v, T const& n) { return root(v, n).id(); }

template<class T>
void set_root(std::vector<T>& v, int g) {
#ifdef _OPENMP
  while(true) {
    int r = root_index(v, g);
    if (r == g) {
      return;
    } else {
      int w = v[r].lock_root();
      if (w != 0) {
        v[g].set_weight(w);
        v[r].set_parent(g); // release lock
        return;
      }
    }
  }
#else
  int r = root_index(v, g);
  if (r != g) {
    v[g] = v[r];
    v[r].set_parent(g);
  }
#endif
}

template<class T>
inline void update_link(std::vector<T>& v, int g, int r) {
  while (g != r) {
    int p = v[g].parent();
    v[g].set_parent(r);
    g = p;
  }
}

// WARNING: this is not thread-safe
template<class T>
inline int unify_compress(std::vector<T>& v, int g0, int g1) {
  using std::swap;
  int r0 = root_index(v, g0);
  int r1 = root_index(v, g1);
  if (r0 != r1) {
#if !defined(CLUSTER_USE_DETERMINISTIC_UNIFY)
    if (v[r0].weight() < v[r1].weight()) swap(r0, r1);
#else
    if (r0 > r1) swap(r0, r1);
#endif
    v[r0].set_weight(v[r0].weight() + v[r1].weight());
    v[r1].set_parent(r0);
  }
  update_link(v, g0, r0);
  update_link(v, g1, r0);
  return r0; // return (new) root node
}

template<class T>
inline int unify_pathhalving(std::vector<T>& v, int g0, int g1) {
  using std::swap;
  int r0 = root_index_ph(v, g0);
  int r1 = root_index_ph(v, g1);
#ifdef _OPENMP
  int w0 = 0;
  int w1 = 0;
  while (true) {
    if (r1 == r0) return r0; // g0 and g1 belong to the same cluster
    w0 = v[r0].lock_root();
    w1 = v[r1].lock_root();
    if (w0 != 0 && w1 != 0) break;
    if (w0 != 0) v[r0].set_weight(w0); // release lock
    if (w1 != 0) v[r1].set_weight(w1); // release lock
    r0 = root_index_ph(v, r0);
    r1 = root_index_ph(v, r1);
  }
# if !defined(CLUSTER_USE_DETERMINISTIC_UNIFY)
  if (w0 < w1) swap(r0, r1);
# else
  if (r0 > r1) swap(r0, r1);
# endif
  v[r0].set_weight(w0 + w1); // release lock
  v[r1].set_parent(r0); // release lock
#else
  if (r0 != r1) {
# if !defined(CLUSTER_USE_DETERMINISTIC_UNIFY)
    if (v[r0].weight() < v[r1].weight()) swap(r0, r1);
# else
    if (r0 > r1) swap(r0, r1);
# endif
    v[r0].set_weight(v[r0].weight() + v[r1].weight());
    v[r1].set_parent(r0);
  }
#endif
  return r0; // return (new) root node
}

template<class T>
inline int unify(std::vector<T>& v, int g0, int g1) {
  return unify_pathhalving(v, g0, g1);
}

template<class T>
inline void output(std::vector<T> const& v, std::ostream& os = std::cout) {
  for (int i = 0; i < v.size(); ++i) {
    os << "node " << i << ": ";
    int g = i;
    if (v[g].is_root()) {
      os << "root (id = " << v[g].id() << ")" << std::endl;
    } else {
      while (true) {
        if (v[g].is_root()) {
          os << g << std::endl;
          break;
        } else {
          os << g << " -> ";
          g = v[g].parent();
        }
      }
    }
  }
}

template<typename T>
int count_root(std::vector<T>& v, int start, int n) {
  int nc = 0;
  for (int i = start; i < start + n; ++i)
    if (v[i].is_root()) ++nc;
  return nc;
}

template<typename T>
int count_root_p(std::vector<T>& v, int start, int n) {
  int nc = 0;
  #pragma omp for schedule(static) nowait
  for (int i = start; i < start + n; ++i)
    if (v[i].is_root()) ++nc;
  return nc;
}

template<typename T>
int set_id(std::vector<T>& v, int start, int n, int nc) {
  for (int i = start; i < start + n; ++i)
    if (v[i].is_root()) v[i].set_id(nc++);
  return nc;
}

template<typename T>
int set_id_p(std::vector<T>& v, int start, int n, int nc) {
  #pragma omp for schedule(static) nowait
  for (int i = start; i < start + n; ++i)
    if (v[i].is_root()) v[i].set_id(nc++);
  return nc;
}

template<typename T>
void copy_id(std::vector<T>& v, int start, int n) {
  for (int i = start; i < start + n; ++i)
    v[i].set_id(cluster_id(v, i));
}

template<typename T>
void copy_id_p(std::vector<T>& v, int start, int n) {
  #pragma omp for schedule(static) nowait
  for (int i = start; i < start + n; ++i)
    v[i].set_id(cluster_id(v, i));
}

template<typename T>
inline void pack_tree(std::vector<T>& v, int n) {
#ifdef _OPENMP
  int g, w; // workaround for FCC OpenMP bug? -- ST 2010-11-22
  #pragma omp parallel for schedule(static) private(g, w)
  for (int i = 0; i < n; ++i) {
    if (!v[i].is_root()) {
      g = v[i].parent();
      while (true) {
        if (g < n) {
          // encounter node with index < n
          v[i].set_parent(g);
          break;
        } else if (v[g].is_root()) {
          // found root with index >= n
          w = v[g].lock_root();
          if (w != 0) {
            v[i].set_weight(w);
            v[g].set_parent(i); // release lock
            break;
          }
        } else {
          g = v[g].parent();
        }
      }
    }
  }
#else
  for (int i = 0; i < n; ++i) {
    if (!v[i].is_root()) {
      int g = v[i].parent();
      while (true) {
        if (g < n) {
          // encounter node with index < n
          v[i].set_parent(g);
          break;
        } else if (v[g].is_root()) {
          // found root with index >= n
          v[i] = v[g];
          v[g].set_parent(i);
          break;
        } else {
          g = v[g].parent();
        }
      }
    }
  }
#endif
}

// pack tree so that nodes with id [0...n) and [m...) come upper
template<typename T>
inline void pack_tree(std::vector<T>& v, int n, int m) {
#ifdef _OPENMP
  int g, w; // workaround for FCC OpenMP bug? -- ST 2010-11-22
  #pragma omp parallel for schedule(static) private(g, w)
  for (int i = 0; i < n; ++i) {
    if (!v[i].is_root()) {
      g = v[i].parent();
      while (true) {
        if (g < n || g >= m) {
          // encounter node with index < n or >= m
          v[i].set_parent(g);
          break;
        } else if (v[g].is_root()) {
          // found root with index >= n and < m
          w = v[g].lock_root();
          if (w != 0) {
            v[i].set_weight(w);
            v[g].set_parent(i); // release lock
            break;
          }
        } else {
          g = v[g].parent();
        }
      }
    }
  }
  #pragma omp parallel for schedule(static) private(g, w)
  for (int i = m; i < m + n; ++i) {
    if (!v[i].is_root()) {
      g = v[i].parent();
      while (true) {
        if (g < n || g >= m) {
          // encounter node with index < n or >= m
          v[i].set_parent(g);
          break;
        } else if (v[g].is_root()) {
          // found root with index >= n and < m
          w = v[g].lock_root();
          if (w != 0) {
            v[i].set_weight(w);
            v[g].set_parent(i); // release lock
            break;
          }
        } else {
          g = v[g].parent();
        }
      }
    }
  }
#else
  for (int i = 0; i < n; ++i) {
    if (!v[i].is_root()) {
      int g = v[i].parent();
      while (true) {
        if (g < n || g >= m) {
          // encounter node with index < n or >= m
          v[i].set_parent(g);
          break;
        } else if (v[g].is_root()) {
          // found root with index >= n and < m
          v[i] = v[g];
          v[g].set_parent(i);
          break;
        } else {
          g = v[g].parent();
        }
      }
    }
  }
  for (int i = m; i < m + n; ++i) {
    if (!v[i].is_root()) {
      int g = v[i].parent();
      while (true) {
        if (g < n || g >= m) {
          // encounter node with index < n or >= m
          v[i].set_parent(g);
          break;
        } else if (v[g].is_root()) {
          // found root with index >= n and < m
          v[i] = v[g];
          v[g].set_parent(i);
          break;
        } else {
          g = v[g].parent();
        }
      }
    }
  }
#endif
}

} // end namespace union_find
} // end namespace cluster
