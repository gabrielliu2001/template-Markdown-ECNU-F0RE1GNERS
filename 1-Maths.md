# Maths

## Linear Sieve

```cpp
void init() {
  phi[1] = 1;
  for (int i = 2; i < MAXN; ++i) {
    if (!vis[i]) {
      pri[cnt++] = i;
    }
    for (int j = 0; j < cnt; ++j) {
      if (1ll * i * pri[j] >= MAXN) break;
      vis[i * pri[j]] = 1;
      if (i % pri[j] == 0) break;
    }
  }
}
```

## Euler Phi Function Sieve

```cpp
void pre() {
  memset(is_prime, 1, sizeof(is_prime));
  int cnt = 0;
  is_prime[1] = 0;
  phi[1] = 1;
  for (int i = 2; i <= 5000000; i++) {
    if (is_prime[i]) {
      prime[++cnt] = i;
      phi[i] = i - 1;
    }
    for (int j = 1; j <= cnt && i * prime[j] <= 5000000; j++) {
      is_prime[i * prime[j]] = 0;
      if (i % prime[j])
        phi[i * prime[j]] = phi[i] * phi[prime[j]];
      else {
        phi[i * prime[j]] = phi[i] * prime[j];
        break;
      }
    }
  }
}
```

## Euler Phi Function

```cpp
auto calc = [&] (int n) {
  int res = n;
  for (int p = 2; p * p <= n; ++p) {
    if (n % p == 0) {
      while (n % p == 0) n /= p;
      res -= res / p;
    }
  }
  if (n > 1) res -= res / n;
  return res;
};
```

## Linear Inverse

```cpp
vl fac, inv, numinv; // ncr and fac

inline ll ncr(int n, int r){
    if (n < 0 || r < 0 || r > n)
        return 0;
    return fac[n] * inv[r] % mod * inv[n - r] % mod;
}
inline void calfacinv(int n){
    fac.reserve(n + 1);
    fac[0] = fac[1] = 1;
    for (int i = 2; i <= n; i++){
        fac[i] = fac[i - 1] * i % mod;
    }
    numinv.reserve(n + 1);
    numinv[0] = numinv[1] = 1;
    for (int i = 2; i <= n; i++){
        numinv[i] = numinv[mod % i] * (mod - mod / i) % mod;
    }
    inv.reserve(n + 1);
    inv[0] = inv[1] = 1;
    for (int i = 2; i <= n; i++){
        inv[i] = numinv[i] * inv[i - 1] % mod;
    }
    return;
}
```

### FFT

+ n 需补成 2 的幂 （n 必须超过 a 和 b 的最高指数之和）

```cpp
typedef double LD;
const LD PI = acos(-1);
struct C {
    LD r, i;
    C(LD r = 0, LD i = 0): r(r), i(i) {}
};
C operator + (const C& a, const C& b) {
    return C(a.r + b.r, a.i + b.i);
}
C operator - (const C& a, const C& b) {
    return C(a.r - b.r, a.i - b.i);
}
C operator * (const C& a, const C& b) {
    return C(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r);
}

void FFT(C x[], int n, int p) {
    for (int i = 0, t = 0; i < n; ++i) {
        if (i > t) swap(x[i], x[t]);
        for (int j = n >> 1; (t ^= j) < j; j >>= 1);
    }
    for (int h = 2; h <= n; h <<= 1) {
        C wn(cos(p * 2 * PI / h), sin(p * 2 * PI / h));
        for (int i = 0; i < n; i += h) {
            C w(1, 0), u;
            for (int j = i, k = h >> 1; j < i + k; ++j) {
                u = x[j + k] * w;
                x[j + k] = x[j] - u;
                x[j] = x[j] + u;
                w = w * wn;
            }
        }
    }
    if (p == -1)
        FOR (i, 0, n)
            x[i].r /= n;
}

void conv(C a[], C b[], int n) {
    FFT(a, n, 1);
    FFT(b, n, 1);
    FOR (i, 0, n)
        a[i] = a[i] * b[i];
    FFT(a, n, -1);
}
```

### FWT

+ $C_k=\sum_{i \oplus j=k} A_i B_j$
+ FWT 完后需要先模一遍

```cpp
template<typename T>
void fwt(LL a[], int n, T f) {
    for (int d = 1; d < n; d *= 2)
        for (int i = 0, t = d * 2; i < n; i += t)
            FOR (j, 0, d)
                f(a[i + j], a[i + j + d]);
}

void AND(LL& a, LL& b) { a += b; }
void OR(LL& a, LL& b) { b += a; }
void XOR (LL& a, LL& b) {
    LL x = a, y = b;
    a = (x + y) % MOD;
    b = (x - y + MOD) % MOD;
}
void rAND(LL& a, LL& b) { a -= b; }
void rOR(LL& a, LL& b) { b -= a; }
void rXOR(LL& a, LL& b) {
    static LL INV2 = (MOD + 1) / 2;
    LL x = a, y = b;
    a = (x + y) * INV2 % MOD;
    b = (x - y + MOD) * INV2 % MOD;
}
```

+ FWT 子集卷积

```text
a[popcount(x)][x] = A[x]
b[popcount(x)][x] = B[x]
fwt(a[i]) fwt(b[i])
c[i + j][x] += a[i][x] * b[j][x]
rfwt(c[i])
ans[x] = c[popcount(x)][x]
```

## simpson 自适应积分

```cpp
LD simpson(LD l, LD r) {
    LD c = (l + r) / 2;
    return (f(l) + 4 * f(c) + f(r)) * (r - l) / 6;
}

LD asr(LD l, LD r, LD eps, LD S) {
    LD m = (l + r) / 2;
    LD L = simpson(l, m), R = simpson(m, r);
    if (fabs(L + R - S) < 15 * eps) return L + R + (L + R - S) / 15;
    return asr(l, m, eps / 2, L) + asr(m, r, eps / 2, R);
}

LD asr(LD l, LD r, LD eps) { return asr(l, r, eps, simpson(l, r)); }
```

## 公式

### 一些数论公式

- 当 $x\geq\phi(p)$ 时有 $a^x\equiv a^{x \; mod \; \phi(p) + \phi(p)}\pmod p$
- $\mu^2(n)=\sum_{d^2|n} \mu(d)$
- $\sum_{d|n} \varphi(d)=n$
- $\sum_{d|n} 2^{\omega(d)}=\sigma_0(n^2)$，其中 $\omega$ 是不同素因子个数
- $\sum_{d|n} \mu^2(d)=2^{\omega(d)}$

### 一些数论函数求和的例子

+ $\sum_{i=1}^n i[gcd(i, n)=1] = \frac {n \varphi(n) + [n=1]}{2}$
+ $\sum_{i=1}^n \sum_{j=1}^m [gcd(i,j)=x]=\sum_d \mu(d) \lfloor \frac n {dx} \rfloor  \lfloor \frac m {dx} \rfloor$
+ $\sum_{i=1}^n \sum_{j=1}^m gcd(i, j) = \sum_{i=1}^n \sum_{j=1}^m \sum_{d|gcd(i,j)} \varphi(d) = \sum_{d} \varphi(d) \lfloor \frac nd \rfloor \lfloor \frac md \rfloor$
+ $S(n)=\sum_{i=1}^n \mu(i)=1-\sum_{i=1}^n \sum_{d|i,d < i}\mu(d) \overset{t=\frac id}{=} 1-\sum_{t=2}^nS(\lfloor \frac nt \rfloor)$
  + 利用 $[n=1] = \sum_{d|n} \mu(d)$
+ $S(n)=\sum_{i=1}^n \varphi(i)=\sum_{i=1}^n i-\sum_{i=1}^n \sum_{d|i,d<i} \varphi(i)\overset{t=\frac id}{=} \frac {i(i+1)}{2} - \sum_{t=2}^n S(\frac n t)$
  + 利用 $n = \sum_{d|n} \varphi(d)$
+ $\sum_{i=1}^n \mu^2(i) = \sum_{i=1}^n \sum_{d^2|n} \mu(d)=\sum_{d=1}^{\lfloor \sqrt n \rfloor}\mu(d) \lfloor \frac n {d^2} \rfloor$ 
+ $\sum_{i=1}^n \sum_{j=1}^n gcd^2(i, j)= \sum_{d} d^2 \sum_{t} \mu(t) \lfloor \frac n{dt} \rfloor ^2 \\
  \overset{x=dt}{=} \sum_{x} \lfloor \frac nx \rfloor ^ 2 \sum_{d|x} d^2 \mu(\frac xd)$
+ $\sum_{i=1}^n \varphi(i)=\frac 12 \sum_{i=1}^n \sum_{j=1}^n [i \perp j] - 1=\frac 12 \sum_{i=1}^n \mu(i) \cdot\lfloor \frac n i \rfloor ^2-1$

### 斐波那契数列性质

- $F_{a+b}=F_{a-1} \cdot F_b+F_a \cdot F_{b+1}$
- $F_1+F_3+\dots +F_{2n-1} = F_{2n},F_2 + F_4 + \dots + F_{2n} = F_{2n + 1} - 1$
- $\sum_{i=1}^n F_i = F_{n+2} - 1$
- $\sum_{i=1}^n F_i^2 = F_n \cdot F_{n+1}$
- $F_n^2=(-1)^{n-1} + F_{n-1} \cdot F_{n+1}$
- $gcd(F_a, F_b)=F_{gcd(a, b)}$
- 模 $n$ 周期（皮萨诺周期）
  - $\pi(p^k) = p^{k-1} \pi(p)$
  - $\pi(nm) = lcm(\pi(n), \pi(m)), \forall n \perp m$
  - $\pi(2)=3, \pi(5)=20$
  - $\forall p \equiv \pm 1\pmod {10}, \pi(p)|p-1$
  - $\forall p \equiv \pm 2\pmod {5}, \pi(p)|2p+2$

### 常见生成函数

+ $(1+ax)^n=\sum_{k=0}^n \binom {n}{k} a^kx^k$
+ $\dfrac{1-x^{r+1}}{1-x}=\sum_{k=0}^nx^k$
+ $\dfrac1{1-ax}=\sum_{k=0}^{\infty}a^kx^k$
+ $\dfrac 1{(1-x)^2}=\sum_{k=0}^{\infty}(k+1)x^k$
+ $\dfrac1{(1-x)^n}=\sum_{k=0}^{\infty} \binom{n+k-1}{k}x^k$
+ $e^x=\sum_{k=0}^{\infty}\dfrac{x^k}{k!}$
+ $\ln(1+x)=\sum_{k=0}^{\infty}\dfrac{(-1)^{k+1}}{k}x^k$

### 佩尔方程

若一个丢番图方程具有以下的形式：$x^2 - ny^2= 1$。且 $n$ 为正整数，则称此二元二次不定方程为**佩尔方程**。

若 $n$ 是完全平方数，则这个方程式只有平凡解 $(\pm 1,0)$（实际上对任意的 $n$，$(\pm 1,0)$ 都是解）。对于其余情况，拉格朗日证明了佩尔方程总有非平凡解。而这些解可由 $\sqrt{n}$ 的连分数求出。

$x = [a_0; a_1, a_2, a_3]=x = a_0 + \cfrac{1}{a_1 + \cfrac{1}{a_2 + \cfrac{1}{a_3 + \cfrac{1}{\ddots\,}}}}$

设 $\tfrac{p_i}{q_i}$ 是 $\sqrt{n}$ 的连分数表示：$[a_{0}; a_{1}, a_{2}, a_{3}, \,\ldots ]$ 的渐近分数列，由连分数理论知存在 $i$ 使得 $(p_i,q_i)$ 为佩尔方程的解。取其中最小的 $i$，将对应的 $(p_i,q_i)$ 称为佩尔方程的基本解，或最小解，记作 $(x_1,y_1)$，则所有的解 $(x_i,y_i)$ 可表示成如下形式：$x_{i}+y_{i}{\sqrt  n}=(x_{1}+y_{1}{\sqrt  n})^{i}$。或者由以下的递回关系式得到：

$\displaystyle x_{i+1} = x_1 x_i + n y_1 y_i$, $\displaystyle y_{{i+1}}=x_{1}y_{i}+y_{1}x_{i}$。

**但是：**佩尔方程千万不要去推（虽然推起来很有趣，但结果不一定好看，会是两个式子）。记住佩尔方程结果的形式通常是 $a_n=ka_{n−1}−a_{n−2}$（$a_{n−2}$ 前的系数通常是 $−1$）。暴力 / 凑出两个基础解之后加上一个 $0$，容易解出 $k$ 并验证。

### Burnside & Polya

+ $|X/G|={\frac  {1}{|G|}}\sum _{{g\in G}}|X^{g}|$

注：$X^g$ 是 $g$ 下的不动点数量，也就是说有多少种东西用 $g$ 作用之后可以保持不变。

+ $|Y^X/G| = \frac{1}{|G|}\sum_{g \in G} m^{c(g)}$

注：用 $m$ 种颜色染色，然后对于某一种置换 $g$，有 $c(g)$ 个置换环，为了保证置换后颜色仍然相同，每个置换环必须染成同色。

### 皮克定理

$2S = 2a+b-2$

+ $S$ 多边形面积
+ $a$ 多边形内部点数
+ $b$ 多边形边上点数

### 莫比乌斯反演

+ $g(n) = \sum_{d|n} f(d) \Leftrightarrow f(n) = \sum_{d|n} \mu (d) g( \frac{n}{d})$
+ $f(n)=\sum_{n|d}g(d) \Leftrightarrow g(n)=\sum_{n|d} \mu(\frac{d}{n}) f(d)$

### 低阶等幂求和

+ $\sum_{i=1}^{n} i^{1} = \frac{n(n+1)}{2} = \frac{1}{2}n^2 +\frac{1}{2} n$
+ $\sum_{i=1}^{n} i^{2} = \frac{n(n+1)(2n+1)}{6} = \frac{1}{3}n^3 + \frac{1}{2}n^2 + \frac{1}{6}n$
+ $\sum_{i=1}^{n} i^{3} = \left[\frac{n(n+1)}{2}\right]^{2} = \frac{1}{4}n^4 + \frac{1}{2}n^3 + \frac{1}{4}n^2$
+ $\sum_{i=1}^{n} i^{4} = \frac{n(n+1)(2n+1)(3n^2+3n-1)}{30} = \frac{1}{5}n^5 + \frac{1}{2}n^4 + \frac{1}{3}n^3 - \frac{1}{30}n$
+ $\sum_{i=1}^{n} i^{5} = \frac{n^{2}(n+1)^{2}(2n^2+2n-1)}{12} = \frac{1}{6}n^6 + \frac{1}{2}n^5 + \frac{5}{12}n^4 - \frac{1}{12}n^2$

### 一些组合公式

+ 错排公式：$D_1=0,D_2=1,D_n=(n-1)(D_{n-1} + D_{n-2})=n!(\frac 1{2!}-\frac 1{3!}+\dots + (-1)^n\frac 1{n!})=\lfloor \frac{n!}e + 0.5 \rfloor$
+ 卡塔兰数（$n$ 对括号合法方案数，$n$ 个结点二叉树个数，$n\times n$ 方格中对角线下方的单调路径数，凸 $n+2$ 边形的三角形划分数，$n$ 个元素的合法出栈序列数）：$C_n=\frac 1{n+1}\binom {2n}n=\frac{(2n)!}{(n+1)!n!}$

## 二次剩余

URAL 1132

```cpp
LL q1, q2, w;
struct P { // x + y * sqrt(w)
    LL x, y;
};

P pmul(const P& a, const P& b, LL p) {
    P res;
    res.x = (a.x * b.x + a.y * b.y % p * w) % p;
    res.y = (a.x * b.y + a.y * b.x) % p;
    return res;
}

P bin(P x, LL n, LL MOD) {
    P ret = {1, 0};
    for (; n; n >>= 1, x = pmul(x, x, MOD))
        if (n & 1) ret = pmul(ret, x, MOD);
    return ret;
}
LL Legendre(LL a, LL p) { return bin(a, (p - 1) >> 1, p); }

LL equation_solve(LL b, LL p) {
    if (p == 2) return 1;
    if ((Legendre(b, p) + 1) % p == 0)
        return -1;
    LL a;
    while (true) {
        a = rand() % p;
        w = ((a * a - b) % p + p) % p;
        if ((Legendre(w, p) + 1) % p == 0)
            break;
    }
    return bin({a, 1}, (p + 1) >> 1, p).x;
}

int main() {
    int T; cin >> T;
    while (T--) {
        LL a, p; cin >> a >> p;
        a = a % p;
        LL x = equation_solve(a, p);
        if (x == -1) {
            puts("No root");
        } else {
            LL y = p - x;
            if (x == y) cout << x << endl;
            else cout << min(x, y) << " " << max(x, y) << endl;
        }
    }
}
```

## 伯努利数和等幂求和

* 预处理逆元
* 预处理组合数
* $\sum_{i=0}^n i^k = \frac{1}{k+1} \sum_{i=0}^k \binom{k+1}{i} B_{k+1-i} (n+1)^i$.
* 也可以 $\sum_{i=0}^n i^k = \frac{1}{k+1} \sum_{i=0}^k \binom{k+1}{i} B^+_{k+1-i} n^i$。区别在于 $B^+_1 =1/2$。(心态崩了)

```cpp
namespace Bernoulli {
    const int M = 100;
    LL inv[M] = {-1, 1};
    void inv_init(LL n, LL p) {
        FOR (i, 2, n)
            inv[i] = (p - p / i) * inv[p % i] % p;
    }

    LL C[M][M];
    void init_C(int n) {
        FOR (i, 0, n) {
            C[i][0] = C[i][i] = 1;
            FOR (j, 1, i)
                C[i][j] = (C[i - 1][j] + C[i - 1][j - 1]) % MOD;
        }
    }

    LL B[M] = {1};
    void init() {
        inv_init(M, MOD);
        init_C(M);
        FOR (i, 1, M - 1) {
            LL& s = B[i] = 0;
            FOR (j, 0, i)
                s += C[i + 1][j] * B[j] % MOD;
            s = (s % MOD * -inv[i + 1] % MOD + MOD) % MOD;
        }
    }

    LL p[M] = {1};
    LL go(LL n, LL k) {
        n %= MOD;
        if (k == 0) return n;
        FOR (i, 1, k + 2)
            p[i] = p[i - 1] * (n + 1) % MOD;
        LL ret = 0;
        FOR (i, 1, k + 2)
            ret += C[k + 1][i] * B[k + 1 - i] % MOD * p[i] % MOD;
        ret = ret % MOD * inv[k + 1] % MOD;
        return ret;
    }
}
```

## 离散对数

### BSGS

+ 模数为素数

```cpp
LL BSGS(LL a, LL b, LL p) { // a^x = b (mod p)
    a %= p;
    if (!a && !b) return 1;
    if (!a) return -1;
    static map<LL, LL> mp; mp.clear();
    LL m = sqrt(p + 1.5);
    LL v = 1;
    FOR (i, 1, m + 1) {
        v = v * a % p;
        mp[v * b % p] = i;
    }
    LL vv = v;
    FOR (i, 1, m + 1) {
        auto it = mp.find(vv);
        if (it != mp.end()) return i * m - it->second;
        vv = vv * v % p;
    }
    return -1;
}
```

### exBSGS

+ 模数可以非素数

```cpp
LL exBSGS(LL a, LL b, LL p) { // a^x = b (mod p)
    a %= p; b %= p;
    if (a == 0) return b > 1 ? -1 : b == 0 && p != 1;
    LL c = 0, q = 1;
    while (1) {
        LL g = __gcd(a, p);
        if (g == 1) break;
        if (b == 1) return c;
        if (b % g) return -1;
        ++c; b /= g; p /= g; q = a / g * q % p;
    }
    static map<LL, LL> mp; mp.clear();
    LL m = sqrt(p + 1.5);
    LL v = 1;
    FOR (i, 1, m + 1) {
        v = v * a % p;
        mp[v * b % p] = i;
    }
    FOR (i, 1, m + 1) {
        q = q * v % p;
        auto it = mp.find(q);
        if (it != mp.end()) return i * m - it->second + c;
    }
    return -1;
}
```

## 数论分块

$f(i) = \lfloor \frac{n}{i} \rfloor=v$ 时 $i$ 的取值范围是 $[l,r]$。

```cpp
for (LL l = 1, v, r; l <= N; l = r + 1) {
    v = N / l; r = N / v;
}
```

## 博弈

+ Nim 游戏：每轮从若干堆石子中的一堆取走若干颗。先手必胜条件为石子数量异或和非零。
+ 阶梯 Nim 游戏：可以选择阶梯上某一堆中的若干颗向下推动一级，直到全部推下去。先手必胜条件是奇数阶梯的异或和非零（对于偶数阶梯的操作可以模仿）。
+ Anti-SG：无法操作者胜。先手必胜的条件是：
  + SG 不为 0 且某个单一游戏的 SG 大于 1 。
  + SG 为 0 且没有单一游戏的 SG 大于 1。
+ Every-SG：对所有单一游戏都要操作。先手必胜的条件是单一游戏中的最大 step 为奇数。
  + 对于终止状态 step 为 0
  + 对于 SG 为 0 的状态，step 是最大后继 step +1
  + 对于 SG 非 0 的状态，step 是最小后继 step +1
+ 树上删边：叶子 SG 为 0，非叶子结点为所有子结点的 SG 值加 1 后的异或和。

尝试：

+ 打表找规律
+ 寻找一类必胜态（如对称局面）
+ 直接博弈 dp