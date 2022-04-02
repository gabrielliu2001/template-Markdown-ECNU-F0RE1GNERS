# Data Structure

## Xor MST

```cpp
struct Trie{
	int son[2][200000*30+10],tot;
	void Insert(int a){
		int now=0,id;
		for(int i=30;i>=0;i--){
	    	id=(a>>i)&1;
	    	if(!son[id][now])son[id][now]=++tot;
	    	now=son[id][now];
		}
	}
	int Find(int r1,int r2,int b){
		if(b<0) return 0;
		int a1=-1,a2=-1;
		if(son[0][r1]&&son[0][r2]) a1=Find(son[0][r1],son[0][r2],b-1);
		if(son[1][r1]&&son[1][r2]) a2=Find(son[1][r1],son[1][r2],b-1);
		if(~a1&&~a2) return min(a1,a2);
		if(~a1) return a1;if(~a2) return a2;
		if(son[1][r1]&&son[0][r2]) a1=Find(son[1][r1],son[0][r2],b-1)+(1<<b);
		if(son[0][r1]&&son[1][r2]) a2=Find(son[0][r1],son[1][r2],b-1)+(1<<b);
		if(~a1&&~a2) return min(a1,a2);
		if(~a1) return a1;if(~a2) return a2;
	}
}T;
long long ans;
void dfs(int a,int b){
	if(b<0) return;
	if(T.son[0][a]&&T.son[1][a]) ans+=1ll*T.Find(T.son[0][a],T.son[1][a],b-1)+(1ll<<b);
	if(T.son[0][a]) dfs(T.son[0][a],b-1);
	if(T.son[1][a]) dfs(T.son[1][a],b-1);
}
int n,v;
int main() {
	n=read();
	for(int i=1;i<=n;i++)T.Insert(read());
	dfs(0,30);
	printf("%I64d\n",ans);
}
```

## First element >= x and index >= l

```cpp
#define maxn 200010
#define ll long long
#define lowbit(i) ((i) & (-i))

int n, m, c[maxn], w[maxn];

ll Bit[maxn];
void add(int i, int v) { while (i <= n) Bit[i] += v, i += lowbit(i); }

ll get_sum(int i) {
    ll s = 0;
    while (i) s += Bit[i], i -= lowbit(i);
    return s; 
}

int pre[maxn];
set<int> S[maxn];

#define lc i << 1
#define rc i << 1 | 1
int T[maxn * 4];
inline void maintain(int i) { T[i] = max(T[lc], T[rc]); } 

void build(int i, int l, int r) {
    if (l == r) return T[i] = pre[l], void();
    int m = l + r >> 1;
    build(lc, l, m); build(rc, m + 1, r);
    maintain(i); 
}

void update(int i, int l, int r, int k, int v) {
    if (l == r) return T[i] = v, void();
    int m = l + r >> 1;
    if (k <= m) update(lc, l, m, k, v);
    else update(rc, m + 1, r, k, v);
    maintain(i);
}

int query(int i, int l, int r, int L, int R, int k) {
    if (l > R || r < L || T[i] < k) return 0;
    if (l == r) return T[i] >= k ? l : 0;
    int m = l + r >> 1, v = query(lc, l, m, L, R, k);
    if (v) return v;
    else return query(rc, m + 1, r, L, R, k); 
}

inline void solve_1() {
    int x, y, z; cin >> x >> y >> z;
    add(x, z - w[x]); w[x] = z;
    auto l = S[c[x]].lower_bound(x), r = S[c[x]].upper_bound(x); --l;
    if (*r != n + 1) pre[*r] = *l, update(1, 1, n, *r, *l);
    S[c[x]].erase(x); c[x] = y; S[c[x]].insert(x); 
    l = S[c[x]].lower_bound(x), r = S[c[x]].upper_bound(x); --l;
    if (*r != n + 1) pre[*r] = x, update(1, 1, n, *r, x);
    pre[x] = *l; update(1, 1, n, x, *l);
}

int tmp[maxn];
inline void solve_2() {
    int x, y; cin >> x >> y;
    vector<int> vec; ll ans = 0; int p = x; 
    while (p <= n) {
        int t = query(1, 1, n, p, n, x);
        if (!t) { ans += get_sum(n) - get_sum(p - 1); break; }
        ans += get_sum(t - 1) - get_sum(p - 1);
        if (!y) break;
        if (!tmp[c[t]]) tmp[c[t]] = w[pre[t]], vec.push_back(c[t]);
        if (tmp[c[t]] < w[t]) ans += w[t] - tmp[c[t]], tmp[c[t]] = w[t];
        p = t + 1; --y;
    } cout << ans << "\n";
    for (auto t : vec) tmp[t] = 0; 
}

int main() {
    cin >> n >> m;
    for (int i = 1; i <= n; ++i) cin >> c[i] >> w[i], S[c[i]].insert(i); 
    for (int i = 1; i <= n; ++i) S[i].insert(0), S[i].insert(n + 1);
    for (int i = 1, last = 0; i <= n; ++i, last = 0)
        for (auto t : S[i]) 
            if (1 <= t && t <= n) pre[t] = last, last = t;
    build(1, 1, n);
    for (int i = 1; i <= n; ++i) add(i, w[i]); 
    for (int i = 1; i <= m; ++i) {
        int opt; cin >> opt;
        if (opt == 1) solve_1();
        else solve_2(); 
    }
}
```

## CDQ

```cpp
#define lowbit(x) ((x)&(-(x)))
const int maxn=100000+10;
int n,m,c[maxn<<1],ans[maxn],cnt;

struct Element{
    int a,b,c,w,f;
}e[maxn],t[maxn];

bool cmp(Element x,Element y){
    if(x.a!=y.a) return x.a<y.a;
    if(x.b!=y.b) return x.b<y.b;
    return x.c<y.c;
}

void update(int x,int y){
    for(;x<=m;x+=lowbit(x)) c[x]+=y;
}
int sum(int x){
    int ans=0;
    for(;x;x-=lowbit(x)) ans+=c[x];
    return ans;
}

void CDQ(int l,int r){
    int mid=(l+r)>>1;
    if(l==r) return ;
    CDQ(l,mid);CDQ(mid+1,r);
    int p=l,q=mid+1,tot=l;
    while(p<=mid&&q<=r){
        if(e[p].b<=e[q].b) update(e[p].c,e[p].w),t[tot++]=e[p++];
        else e[q].f+=sum(e[q].c),t[tot++]=e[q++];
    }
    while(p<=mid) update(e[p].c,e[p].w),t[tot++]=e[p++];
    while(q<=r) e[q].f+=sum(e[q].c),t[tot++]=e[q++];
    for(int i=l;i<=mid;i++) update(e[i].c,-e[i].w);
    for(int i=l;i<=r;i++) e[i]=t[i];
}

int main() {
    n=read();m=read();
    for(int i=1;i<=n;i++)
        e[i].a=read(),e[i].b=read(),e[i].c=read(),e[i].w=1;
    sort(e+1,e+n+1,cmp);
    cnt=1;
    for(int i=2;i<=n;i++){
        if(e[i].a==e[cnt].a&&e[i].b==e[cnt].b&&e[i].c==e[cnt].c) e[cnt].w++;
        else e[++cnt]=e[i];
    }
    CDQ(1,cnt);
    for(int i=1;i<=cnt;i++) ans[e[i].f+e[i].w-1]+=e[i].w;
    for(int i=0;i<n;i++) printf("%d\n",ans[i]);
}
```

## Parallel Binary Search

```cpp
#define lowbit(x) ((x)&(-(x)))
const int maxn=200000+10;
const int inf=1e9;
int n,m,a[maxn],c[maxn],ans[maxn],cnt,tot;

struct Query{
    int l,r,k,id,op;
}q[maxn*3],q1[maxn*3],q2[maxn*3];

void add(int x,int y){
    for(;x<=n;x+=lowbit(x)) c[x]+=y;
}
int sum(int x){
    int ans=0;
    for(;x;x-=lowbit(x)) ans+=c[x];
    return ans;
}

void solve(int l,int r,int L,int R){
    if(L > R) return ;
    if(l == r){
        for(int i=L;i<=R;i++) 
            if(q[i].op==2) ans[q[i].id]=l;
        return ; 
    }
    int mid=(l+r)>>1,cnt1=0,cnt2=0,x;
    for(int i=L;i<=R;i++){
        if(q[i].op==1){
            if(q[i].l <= mid) q1[++cnt1]=q[i],add(q[i].id,q[i].r);
            else q2[++cnt2]=q[i];
        }
        else {
            x=sum(q[i].r)-sum(q[i].l-1);
            if(q[i].k <= x) q1[++cnt1]=q[i];
            else q[i].k-=x,q2[++cnt2]=q[i];
        }
    }
    for(int i=1;i<=cnt1;i++)
        if(q1[i].op==1) add(q1[i].id,-q1[i].r);
    for(int i=1;i<=cnt1;i++) q[L+i-1]=q1[i];
    for(int i=1;i<=cnt2;i++) q[L+i+cnt1-1]=q2[i];
    solve(l,mid,L,L+cnt1-1);
    solve(mid+1,r,L+cnt1,R);
}

int main() {
    n=read(),m=read();
    int l,r,k;char op;
    for(int i=1;i<=n;i++) a[i]=read(),q[++cnt]=(Query){a[i],1,0,i,1};
    for(int i=1;i<=m;i++){
        op=getchar();
        while(!isalpha(op)) op=getchar();
        if(op=='Q') l=read(),r=read(),k=read(),q[++cnt]=(Query){l,r,k,++tot,2};
        else l=read(),r=read(),q[++cnt]=(Query){a[l],-1,0,l,1},q[++cnt]=(Query){a[l]=r,1,0,l,1};
    }
    solve(-inf,inf,1,cnt);
    for(int i=1;i<=tot;i++) printf("%d\n",ans[i]);
    return 0;
}
```

## Segment Tree D & C

```cpp
const int N = 1e5 + 7, M = 2e5 + 7;
int n, m, k, u[M], v[M], f[N<<1], d[N<<1];
struct T {
    int l, r;
    vi e;
} t[N<<2];
stack< pi > s;

void build(int p, int l, int r) {
    t[p].l = l, t[p].r = r;
    if (l == r) return;
    build(ls, l, md), build(rs, md + 1, r);
}

void ins(int p, int l, int r, int x) {
    if (t[p].l >= l && t[p].r <= r) return t[p].e.pb(x), void();
    if (l <= md) ins(ls, l, r, x);
    if (r > md) ins(rs, l, r, x);
}

inline int get(int x) {
    while (x ^ f[x]) x = f[x];
    return x;
}

inline void merge(int x, int y) {
    if (x == y) return;
    if (d[x] > d[y]) swap(x, y);
    s.push(mp(x, d[x] == d[y])), f[x] = y, d[y] += d[x] == d[y];
}

void dfs(int p, int l, int r) {
    bool ok = 1;
    ui o = s.size();
    for (ui i = 0; i < t[p].e.size(); i++) {
        int x = t[p].e[i], u = get(::u[x]), v = get(::v[x]);
        if (u == v) {
            for (int j = l; j <= r; j++) prints("No");
            ok = 0;
            break;
        }
        merge(get(::u[x] + N), v), merge(get(::v[x] + N), u);
    }
    if (ok) {
        if (l == r) prints("Yes");
        else dfs(ls, l, md), dfs(rs, md + 1, r);
    }
    while (s.size() > o) d[f[s.top().fi]] -= s.top().se, f[s.top().fi] = s.top().fi, s.pop();
}

int main() {
    rd(n), rd(m), rd(k), build(1, 1, k);
    for (int i = 1, l, r; i <= m; i++) {
        rd(u[i]), rd(v[i]), rd(l), rd(r);
        if (l ^ r) ins(1, l + 1, r, i);
    }
    for (int i = 1; i <= n; i++) f[i] = i, f[i+N] = i + N;
    dfs(1, 1, k);
    return 0;
}
```

## Mo's Algorithm

```cpp
int main() {
  int n;
  cin >> n;
  vector<int> a(n);
  for (auto& o : a) {
    cin >> o;
    --o;
  }
  int q;
  cin >> q;
  vector queries(q, tuple(0, 0, 0));
  for (int i = 0; i < q; ++i) {
    int l, r;
    cin >> l >> r;
    queries[i] = {l - 1, r, i};
  }
  const int BLOCK_SIZE = int(n / sqrt(q)) > 0 ? n / sqrt(q) : sqrt(n);
  sort(queries.begin(), queries.end(), [&](const auto& x, const auto& y) {
    auto [xl, xr, _i] = x;
    auto [yl, yr, _j] = y;
    if (xl / BLOCK_SIZE != yl / BLOCK_SIZE) {
      return xl / BLOCK_SIZE < yl / BLOCK_SIZE;
    }
    return (xl / BLOCK_SIZE) & 1 ? xr > yr : xr < yr;
  });
  vector<int> ans(q), cnt(n);
  int cur = 0;
  auto add = [&](int i, int v) {
    cnt[a[i]] += v;
    cur += v * ((cnt[a[i]] & 1) == (v < 0));
  };
  // [l, r)
  for (int l = 0, r = 0; auto [ql, qr, i] : queries) {
    while (l > ql) add(--l, 1);
    while (r < qr) add(r++, 1);
    while (l < ql) add(l++, -1);
    while (r > qr) add(--r, -1);
    ans[i] = cur;
  }
  for (auto o : ans) {
    cout << o << '\n';
  }
}
```

## Link/Cut Tree
```cpp
#include<bits/stdc++.h>
#define R register int
#define I inline void
#define G if(++ip==ie)if(fread(ip=buf,1,SZ,stdin))
#define lc c[x][0]
#define rc c[x][1]
using namespace std;
const int SZ=1<<19,N=3e5+9;
char buf[SZ],*ie=buf+SZ,*ip=ie-1;
inline int in(){
	G;while(*ip<'-')G;
	R x=*ip&15;G;
	while(*ip>'-'){x*=10;x+=*ip&15;G;}
	return x;
}
int f[N],c[N][2],v[N],s[N],st[N];
bool r[N];
inline bool nroot(R x){//判断节点是否为一个Splay的根（与普通Splay的区别1）
	return c[f[x]][0]==x||c[f[x]][1]==x;
}//原理很简单，如果连的是轻边，他的父亲的儿子里没有它
I pushup(R x){//上传信息
	s[x]=s[lc]^s[rc]^v[x];
}
I pushr(R x){R t=lc;lc=rc;rc=t;r[x]^=1;}//翻转操作
I pushdown(R x){//判断并释放懒标记
	if(r[x]){
		if(lc)pushr(lc);
		if(rc)pushr(rc);
		r[x]=0;
	}
}
I rotate(R x){//一次旋转
	R y=f[x],z=f[y],k=c[y][1]==x,w=c[x][!k];
	if(nroot(y))c[z][c[z][1]==y]=x;c[x][!k]=y;c[y][k]=w;//额外注意if(nroot(y))语句，此处不判断会引起致命错误（与普通Splay的区别2）
	if(w)f[w]=y;f[y]=x;f[x]=z;
	pushup(y);
}
I splay(R x){//只传了一个参数，因为所有操作的目标都是该Splay的根（与普通Splay的区别3）
	R y=x,z=0;
	st[++z]=y;//st为栈，暂存当前点到根的整条路径，pushdown时一定要从上往下放标记（与普通Splay的区别4）
	while(nroot(y))st[++z]=y=f[y];
	while(z)pushdown(st[z--]);
	while(nroot(x)){
		y=f[x];z=f[y];
		if(nroot(y))
			rotate((c[y][0]==x)^(c[z][0]==y)?x:y);
		rotate(x);
	}
	pushup(x);
}
/*当然了，其实利用函数堆栈也很方便，代替上面的手工栈，就像这样
I pushall(R x){
	if(nroot(x))pushall(f[x]);
	pushdown(x);
}*/
I access(R x){//访问
	for(R y=0;x;x=f[y=x])
		splay(x),rc=y,pushup(x);
}
I makeroot(R x){//换根
	access(x);splay(x);
	pushr(x);
}
int findroot(R x){//找根（在真实的树中的）
	access(x);splay(x);
	while(lc)pushdown(x),x=lc;
	splay(x);
	return x;
}
I split(R x,R y){//提取路径
	makeroot(x);
	access(y);splay(y);
}
I link(R x,R y){//连边
	makeroot(x);
	if(findroot(y)!=x)f[x]=y;
}
I cut(R x,R y){//断边
	makeroot(x);
	if(findroot(y)==x&&f[y]==x&&!c[y][0]){
		f[y]=c[x][1]=0;
		pushup(x);
	}
}
int main() {
	R n=in(),m=in();
	for(R i=1;i<=n;++i)v[i]=in();
	while(m--){
		R type=in(),x=in(),y=in();
		switch(type){
		case 0:split(x,y);printf("%d\n",s[y]);break;
		case 1:link(x,y);break;
		case 2:cut(x,y);break;
		case 3:splay(x);v[x]=y;//先把x转上去再改，不然会影响Splay信息的正确性
		}
	}
	return 0;
}
```

## Filler

```cpp
const int N = 205;
ll n, m, Q;
ll t1[N][N], t2[N][N], t3[N][N], t4[N][N];

void add(ll x, ll y, ll z){
  for(int X = x; X <= n; X += X & -X)
    for(int Y = y; Y <= m; Y += Y & -Y){
      t1[X][Y] += z;
      t2[X][Y] += z * x;
      t3[X][Y] += z * y;
      t4[X][Y] += z * x * y;
    }
}

void range_add(ll xa, ll ya, ll xb, ll yb, ll z){ //(xa, ya) 到 (xb, yb) 的矩形
  add(xa, ya, z);
  add(xa, yb + 1, -z);
  add(xb + 1, ya, -z);
  add(xb + 1, yb + 1, z);
}

ll ask(ll x, ll y){
  ll res = 0;
  for(int i = x; i; i -= i & -i)
  for(int j = y; j; j -= j & -j)
    res += (x + 1) * (y + 1) * t1[i][j]
          - (y + 1) * t2[i][j]
          - (x + 1) * t3[i][j]
          + t4[i][j];
  return res;
}

ll range_ask(ll xa, ll ya, ll xb, ll yb){
  return ask(xb, yb) - ask(xb, ya - 1) - ask(xa - 1, yb) + ask(xa - 1, ya - 1);
}
```

## Two pointers without deletion

```cpp
typedef long long ll;
const int MAXN = 200005;
int T, N;
ll A[MAXN], resl[MAXN], resr;
int main() {
	for (scanf("%d", &T); T; T--) {
		scanf("%d", &N);
		for (int i=1; i<=N; i++) scanf("%lld", &A[i]);
		if (N==1) { puts("1"); continue; }
		for (int i=2; i<=N; i++) A[i-1] -= A[i], A[i-1] = abs(A[i-1]);
		//for (int i=1; i< N; i++) printf("%lld ", A[i]); puts("");
		int l = 1, r = 1, mid = 1, ans = A[1] > 1; // [l, mid], (mid, r];
		resl[1] = A[1]; if (A[1]==1) l = mid+1;
		while (r< N-1) {
			++r, resr = r==mid+1 ? A[r] : gcd(resr, A[r]);
			while (l<=mid && gcd(resl[l], resr)==1) l++;
			if (l> mid) {
				mid = r, l = r+1, resl[l] = A[l-1];
				while (l> 1 && (resl[l-1]=gcd(resl[l], A[l-1]))> 1) l--;
			}
			ans = max(ans, r-l+1);
			//printf("[%d, %d]\n", l, r);
		}
		printf("%d\n", ans+1);
	}
}
```