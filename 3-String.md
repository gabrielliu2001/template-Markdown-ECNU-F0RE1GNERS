# String

## Generalized Suffix Automaton

```cpp
namespace SA {
  constexpr int N = 4e5;
  int ch[N][26], vis[N][26];
  int last, tot, len[N], link[N];

  inline void init() {
    last = 0;
    tot = 1;
    link[0] = -1;
  }

  inline void extend(int c) {
    int u = last;
    if (ch[u][c]) {
      if (len[u] + 1 == len[ch[u][c]]) {
        last = ch[u][c];
        return;
      }
      int w = ch[u][c], clone = tot++;
      memcpy(ch[clone], ch[w], sizeof ch[w]);
      len[clone] = len[u] + 1;
      link[clone] = link[w];
      link[w] = clone;
      for (; u != -1 && ch[u][c] == w; u = link[u]) {
        ch[u][c] = clone;
      }
      last = clone;
      return;
    }
    int v = tot++;
    len[v] = len[u] + 1;
    for (; u != -1 && !ch[u][c]; u = link[u]) {
      ch[u][c] = v;
    }
    if (u == -1) {
      link[v] = 0;
    } else if (len[u] + 1 == len[ch[u][c]]) {
      link[v] = ch[u][c];
    } else {
      int w = ch[u][c], clone = tot++;
      memcpy(ch[clone], ch[w], sizeof ch[w]);
      len[clone] = len[u] + 1;
      link[clone] = link[w];
      link[w] = link[v] = clone;
      for (; u != -1 && ch[u][c] == w; u = link[u]) {
        ch[u][c] = clone;
      }
    }
    last = v;
  }
}

int main() {
  int n;
  cin >> n;
  vector<string> s(n);
  vector<int> ec(n);
  using namespace SA;
  init();
  for (int i = 0; i < n; ++i) {
    cin >> s[i];
    last = 0;
    for (auto& c : s[i]) {
      extend(c - 'a');
    }
    ec[i] = last;
  }
  for (int i = 0; i < n; ++i) {
    int u = 0;
    for (auto c : s[i]) {
      c -= 'a';
      vis[u][c] = 1;
      u = ch[u][c];
    }
  }
  string t;
  cin >> t;
  int u = 0;
  vector<int> match(tot), cnt(tot), q(tot);
  for (auto c : t) {
    c -= 'a';
    for (; u != -1 && !(ch[u][c] && vis[u][c]); u = link[u]);
    if (u == -1) {
      u = 0;
      continue;
    }
    u = ch[u][c];
    ++match[u];
  }
  for (int i = 0; i < tot; ++i) {
    ++cnt[len[i]];
  }
  for (int i = 1; i < tot; ++i) {
    cnt[i] += cnt[i - 1];
  }
  for (int i = 0; i < tot; ++i) {
    q[--cnt[len[i]]] = i;
  }
  for (int i = tot - 1; i > 0; --i) {
    match[link[q[i]]] += match[q[i]];
  }
  for (auto& o : ec) {
    cout << match[o] << '\n';
  }
}
```