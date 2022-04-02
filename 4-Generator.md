# Generator

## makefile

```makefile
% : %.cpp
	g++ -g -std=c++2a $< -o $@ -O2

# % : %.cpp
# 	g++ -Wall -Wconversion -Wfatal-errors -Wshadow -g -std=c++17 -fsanitize=undefined,address $< -o $@ -O2
```

## bash script
```sh
#!/bin/bash

cnt=1
cnt_max=12
wa=0
g++ -std=c++2a gen.cpp -o gen -O2
g++ -std=c++2a all.cpp -o all -O2
g++ -std=c++2a all2.cpp -o all2 -O2
# make gen
# make all
# make all2

# im=Impossible
# im=-1

# while [ $wa -eq 0 ] && [ $cnt -le $cnt_max ]
while [ $wa -eq 0 ]
do
    ./gen ${cnt} > all.in
    # echo Case# $cnt input:
    ./all < all.in > all.out
    # cat "all.in" > "all2.in"
    # cat "all.out" >> "all2.in"
    ./all2 < all.in > all2.out
    # python3 all.py < all.in > all2.out
    diff all.out all2.out
    wa=$?
    [ $wa == 1 ] && s="WA" || s="AC"
    echo Case# $cnt result: $s
    # line=$(head -n 1 all.out)
    # if [ "$line" != "$im" ]
    # then echo Case# $cnt answer: 
        # cat all.out
    # then cp all.in ${cnt}.in
    # fi
    cnt=$(($cnt + 1))
done

# zip -m input.zip *[0-9].in
```

## minified "testlib.h"

```cpp
// ICPC Generator by gabrielliu2001 (Mar 31 2022)
#include "bits/stdc++.h"
using namespace std;

// https://github.com/MikeMirzayanov/testlib
// Adapted from "testlib.h" for stressing solutions in ICPC
struct icpc_gen {
  unsigned long long seed = 3905348978240129619LL;
  unsigned long long mul = 0x5DEECE66DLL;
  unsigned long long add = 0xBLL;
  unsigned long long mask = (1LL << 48) - 1;
  mt19937_64 rng;

  void setSeed(int argc, char *argv[], int fix_seed) {
    if (fix_seed) {
      for (int i = 1; i < argc; ++i) {
        size_t le = strlen(argv[i]);
        for (size_t j = 0; j < le; ++j) {
          seed = seed * mul + (unsigned int) (argv[i][j]) + add;
        }
        seed += mul / add;
      }
    } else {
      seed = chrono::steady_clock::now().time_since_epoch().count();
    }
    rng.seed(seed);
  }

  // Random value in range [l, r] for int, range [l, r) for real
  template<class T>
  T next(T l, T r) {
    assert(l <= r);
    if constexpr(is_floating_point_v<T>) {
      return uniform_real_distribution<T>(l, r)(rng);
    } else {
      return uniform_int_distribution<T>(l, r)(rng);
    }
  }

  // Random permutation of the given size (values are between `f` and `f + n - 1`)
  template<class T, class E = int>
  vector<E> perm(T n, E f = 0) {
    assert(n > 0);
    vector<E> p(n);
    iota(p.begin(), p.end(), f);
    shuffle(p.begin(), p.end(), rng);
    return p;
  }

  // Return `size` unordered (unsorted) distinct numbers between `l` and `r`
  template<class T>
  vector<T> distinct(int size, T l, T r) {
    vector<T> v;
    assert(size > 0);
    assert(l <= r);
    T n = r - l + 1;
    assert(size <= n);
    double ev = 0.0;
    for (int i = 1; i <= size; ++i) {
      ev += double(n) / double(n - i + 1);
    }
    if (ev < double(n)) {
      set<T> s;
      while (int(s.size()) < size) {
        T x = T(next(l, r));
        if (s.emplace(x).second) {
          v.emplace_back(x);
        }
      }
    } else {
      assert(n <= int(1e9));
      vector<T> p(perm(int(n), l));
      v.insert(v.end(), p.begin(), p.begin() + size);
    }
    return v;
  }

  // Return random (unsorted) partition which is a representation of sum as a sum of integers not less than min_part
  template<class T>
  vector<T> partition(int size, T sum, T min_part = 1) {
    assert(size > 0);
    assert(min_part * size <= sum);
    T given_sum = sum, result_sum = 0;
    sum -= min_part * size;
    vector<T> septums(size), result(size);
    vector<T> d = distinct(size - 1, T(1), T(sum + size - 1));
    for (int i = 0; i + 1 < size; ++i) {
      septums[i + 1] = d[i];
    }
    sort(septums.begin(), septums.end());
    for (int i = 0; i + 1 < size; ++i) {
      result[i] = septums[i + 1] - septums[i] - 1;
    }
    result[size - 1] = sum + size - 1 - septums.back();
    for (auto& o : result) {
      o += min_part;
      result_sum += o;
    }
    assert(result_sum == given_sum);
    assert(*min_element(result.begin(), result.end()) >= min_part);
    assert(int(result.size()) == size && result.size() == (size_t) size);
    return result;
  }

  // Random interval [l, r] where l <= r
  template<class T>
  pair<T, T> interval(T l, T r) {
    assert(l <= r);
    T x = next(l - 1, r);
    T y = next(l, r);
    return x == l - 1 ? pair(y, y) : pair(min(x, y), max(x, y));
  }
} rnd;

int main(int argc, char* argv[]) {
  ios_base::sync_with_stdio(0), cin.tie(0);
  rnd.setSeed(argc, argv, 1);
}
```

## Tree Generator

```cpp
// Tree Generator by gabrielliu2001 (Mar 31 2022)
#include "testlib.h"
#include "bits/stdc++.h"
using namespace std;

struct tree_generator {
  // https://cp-algorithms.com/graph/pruefer_code.html#restoring-the-tree-using-the-prufer-code-in-linear-time
  // Random tree by decoding random generated pruefer code
  vector<pair<int, int>> gen_random_tree(int n) {
    if (n == 1) {
      return vector<pair<int, int>>();
    }
    vector<int> code(n - 2);
    for (int& v : code) v = rnd.next(n);

    vector<int> degree(n, 1);
    for (int i : code) degree[i]++;

    int ptr = 0;
    while (degree[ptr] != 1) ptr++;
    int leaf = ptr;

    vector<pair<int, int>> edges;
    for (int v : code) {
      edges.emplace_back(leaf, v);
      if (--degree[v] == 1 && v < ptr) {
        leaf = v;
      } else {
        ptr++;
        while (degree[ptr] != 1) ptr++;
        leaf = ptr;
      }
    }
    edges.emplace_back(leaf, n - 1);
    return edges;
  }

  // i <= chain_sz --> chain, i > chain_sz --> center
  // is_caterpillar --> random center
  // strict_caterpillar --> single node connect to chain node only
  vector<pair<int, int>> gen_chain_related(int n, int chain_sz, bool is_caterpillar = false, bool strict_caterpillar = false) {
    vector<pair<int, int>> edges;
    vector<int> f = rnd.perm(n);
    for (int i = 1; i < n; ++i) {
      int p = i <= chain_sz ? i - 1 : (is_caterpillar ? (strict_caterpillar ? rnd.next(0, chain_sz) : rnd.next(i)) : 0);
      edges.emplace_back(f[p], f[i]);
    }
    return edges;
  }

  // sz = number of children a parent has
  // node_sz = actual number of vertices within a node
  vector<pair<int, int>> gen_binary_tree_related(int n, int sz, int node_sz = 1) {
    vector<pair<int, int>> edges;
    vector<vector<int>> nodes((n - 1) / node_sz + 1);
    vector<int> f = rnd.perm(n);
    for (int i = 0; i < n; ++i) {
      nodes[i / node_sz].emplace_back(i);
    }
    for (auto v : nodes) {
      for (int i = 1; i < v.size(); ++i) {
        edges.emplace_back(f[v[i - 1]], f[v[i]]);
      }
    }
    for (int i = 1; i < nodes.size(); ++i) {
      int u = nodes[(i - 1) / sz].back(), v = nodes[i].front();
      edges.emplace_back(f[u], f[v]);
    }
    return edges;
  }

  // Binary chain
  vector<pair<int, int>> gen_binary_chain(int n) {
    vector<pair<int, int>> edges;
    vector<int> f = rnd.perm(n);
    vector<int> tri;
    for (int acc = 1, sum = 1; sum <= n; ++acc, sum += acc) {
      tri.emplace_back(sum);
    }
    for (int i = 1; i < tri.size(); ++i) {
      for (int j = tri[i - 1]; j < tri[i]; ++j) {
        int p = j - i - (j + 1 == tri[i]);
        edges.emplace_back(f[p], f[j]);
      }
    }
    for (int i = tri.back(); i < n; ++i) {
      int p = rnd.next(i);
      edges.emplace_back(f[p], f[i]);
    }
    return edges;
  }

  // https://oi-wiki.org/contest/problemsetting/#_26
  // type 0: random tree, 1: chain, 2: star, 3: (n / 2) chain + (n / 2) star
  // {4, 5, 6}: caterpillar with {[5, 10], sqrt(n), (n / 2)} single nodes
  // 7: binary tree, 8: binary sqrt(n) chain, 9: sqrtary tree, 10: sqrtray sqrt(n) chain
  // 11: binary chain of depth d --> chain of depth d and binary chain of depth (d - 1)
  vector<pair<int, int>> gen_common_tree(int n, int type = 0, int base = 0) {
    vector<pair<int, int>> edges;
    switch (type) {
      case 0: edges = gen_random_tree(n); break;
      case 1: edges = gen_chain_related(n, n); break;
      case 2: edges = gen_chain_related(n, 0); break;
      case 3: edges = gen_chain_related(n, n / 2); break;
      case 4: edges = gen_chain_related(n, n - rnd.next(5, 10), true); break;
      case 5: edges = gen_chain_related(n, n - int(sqrt(n)), true); break;
      case 6: edges = gen_chain_related(n, n / 2, true); break;
      case 7: edges = gen_binary_tree_related(n, 2); break;
      case 8: edges = gen_binary_tree_related(n, 2, sqrt(n)); break;
      case 9: edges = gen_binary_tree_related(n, sqrt(n)); break;
      case 10: edges = gen_binary_tree_related(n, sqrt(n), sqrt(n)); break;
      case 11: edges = gen_binary_chain(n); break;
      default: edges = gen_random_tree(n);
    }
    for (auto& [u, v] : edges) {
      u += base, v += base;
      if (rnd.next(2)) {
        swap(u, v);
      }
    }
    shuffle(edges.begin(), edges.end());
    return edges;
  }
} tg;

int main(int argc, char* argv[]) {
  registerGen(argc, argv, 1);
  int n = opt<int>("n");
  int type = opt<int>("type");
  int base = opt<int>("base");
  
  vector<pair<int, int>> edges = tg.gen_common_tree(n, type, base);

  println(n);
  for (auto [u, v] : edges) {
    println(u, v);
  }
}
```

## Prime Related Generator

```py
# ICPC Generator (Prime Related) by gabrielliu2001 (Mar 31 2022)

# Maximize number of duplicate prime divisor: power of 2
# Maximize number of prime divisor: 2 * 3 * 5 * 7 * ...
# Highly composite numbers:
def divisor_count(n):
  i = 2
  cnt = 0
  while i ** 2 <= n:
    if n % i == 0:
      cnt += 2
      if n // i == i:
        cnt -= 1
    i += 1
  return cnt

A002182_list, r = [], 0
i = 1
while i <= 10 ** 5:
  d = divisor_count(i)
  if d > r:
    r = d
    A002182_list.append(i)
  i += 1

print(A002182_list)
```