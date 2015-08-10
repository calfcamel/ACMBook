# FFT


---

##复数系

[已经验证(hdu1402 A * B Problem Plus)][1]
```C++
struct complex
{
    double r,i;
    complex(double r = 0,double i = 0):r(r),i(i){};
    complex operator + (const complex &b){return complex(r + b.r, i + b.i);}
    complex operator - (const complex &b){return complex(r - b.r, i - b.i);}
    complex operator * (const complex &b){return complex(r * b.r - i * b.i, r * b.i + i * b.r);}
};

void change(complex y[],int len)
{
    int i,j,k;
    for(i = 1, j = len / 2; i < len - 1; i++)
    {
        if(i < j) swap(y[i], y[j]);
        k = len >> 1;
        while( j >= k)
        {
            j -= k;
            k >>= 1;
        }
        if(j < k)j += k;
    }
}

void fft(complex y[],int len,int on)
{
    change(y,len);
    for(int h = 2;h <= len;h <<= 1)
    {
        complex wn(cos(-on*2*PI/h),sin(-on*2*PI/h));
        for(int j = 0;j < len;j += h)
        {
            complex w(1,0);
            for(int k = j;k < j+h/2;k++)
            {
                complex u = y[k];
                complex t = w*y[k+h/2];
                y[k] = u+t;
                y[k+h/2] = u-t;
                w = w*wn;
            }
        }
    }
    if(on == -1)
        for(int i = 0;i < len;i++)
            y[i].r /= len;
}

//调用
    fft(x1,l,1);
    fft(x2,l,1);
    for(i = 0; i < l; i++)
        x1[i] = x1[i] * x2[i];
    fft(x1,l,-1);
```

---

##大整数系

```C++
// FFT 大整数板子
// 先mk
// 然后conv

#include <cstring>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <stack>
#include <string>
#include <bitset>
using namespace std;
typedef long long LL;
#define clr(a,b) memset(a,b,sizeof(a))
const int MAXN = 100000 + 5;
const int ADD = 50000 + 5;
const LL MOD = 1000000000 + 7;

const double PI = acos(-1.0);


///====================================================
///====================================================
//const LL P = 479LL << 21LL ^ 1LL;
const LL P = 50000000001507329LL;
//const LL P = 998244353LL;
const LL G = 3;
LL g[100];
LL BIT_CNT;
LL N;
LL Mul(LL x, LL y) {
    return ((x*y-(long long)(x/(long double)P*y+1e-8)*P)%P+P)%P;
}
LL pw(LL x, LL y)
{
    LL ret = 1;
    while(y)
    {
        if(y & 1LL) ret = Mul(ret, x);
        y >>= 1LL;
        x = Mul(x,x);
    }
    return ret;
}
inline LL reverse(LL x)
{
    LL ret = 0;
    for(LL i = 0; i < BIT_CNT; i++)
        if(x & (1LL << i)) ret |= 1LL << (BIT_CNT - i - 1LL);
    return ret;
}
inline void mk(LL *x,LL *y)
{
    for(LL i = 0; i < N; i ++)
        x[reverse(i)] = y[i];
}
void FFT(LL *x){
    LL m,i,j,l=1,tk,tx;long long t;
    for(m=-1;++m!=BIT_CNT;l<<=1)for(i=-l<<1;(i+=l<<1)!=N;)for(j=-1,t=1;++j!=l;t=Mul(g[m],t)){
        x[i+j]=((tx=x[i+j])+(tk=Mul(t,x[i+j+l])))%P;
        if((x[i+j+l]=(tx-tk)%P)<0)x[i+j+l]+=P;
    }
}

void conv(LL *x, LL *y) //calc x * y
{
    LL t = P - 1LL;
    for(LL i = 0; i < BIT_CNT; i++)
    {
        t >>= 1;
        g[i] = pw(G,t);
    }
    FFT(x);
    FFT(y);
    for(LL i = 0; i < N; i ++) y[i] = Mul(x[i], y[i]);
    for(LL i = 0; i < N; i ++) x[reverse(i)] = y[ i == 0 ? 0 : N - i ];
    FFT(x);
    LL inv = pw(N,P - 2);
    for(LL i = 0; i < N; i ++) x[i] = Mul(x[i], inv);
}

///====================================================
///====================================================

LL ans[MAXN << 2];

LL a[MAXN];
LL s[MAXN];
LL xx[MAXN << 2];
LL x1[MAXN << 2];
LL x2[MAXN << 2];
int main()
{
    //freopen("h.in","r",stdin);freopen("hh.out","w",stdout);
    int _T;
    scanf("%d",&_T);
    while(_T--)
    {
        LL n; scanf("%lld",&n);
        LL SUM = 0;
        for(LL i = 0; i < n; i++)
        {
            scanf("%lld",&a[i]);
            SUM += a[i];
            s[i] = SUM;
        }
        N = 1;
        N = 1;
        BIT_CNT = 0;
        while(N < SUM + SUM + 2) {N <<= 1LL;BIT_CNT++;}
        //cout << "SUM = " << SUM << " N = " << N << " BIT_CNT = " << BIT_CNT << endl;
        clr(xx,0);
        for(LL i = 0; i < n; i++)
            xx[s[i]] += i + 1;
        //for(LL i = 0; i < N; i++)
        //    cout << xx[i] << " "; cout << endl;
        mk(x1,xx);
        clr(xx,0);
        for(LL i = 0; i < n - 1; i++)
            xx[SUM - s[i]] ++;
        xx[SUM]++;
        //for(LL i = 0; i < N; i++)
        //    cout << xx[i] << " "; cout << endl;
        mk(x2,xx);
        conv(x1,x2);
        for(LL i = SUM; i < SUM + SUM + 1; i++)
            ans[i] = x1[i];
        //cout << "cout 1" << endl;
        //for(LL i = 0; i < N; i++)
        //    cout << ans[i] << " "; cout << endl;
        //cout << "END cout 1" << endl;
        //for(LL i = SUM; i < SUM + SUM + 1; i++)
        //    cout << ans[i] << " "; cout << endl;

        clr(xx,0);
        for(LL i = 0; i < n; i++)
            xx[s[i]] += 1;
        mk(x1,xx);
        clr(xx,0);
        for(LL i = 0; i < n - 1; i++)
            xx[SUM - s[i]] += i + 1;
        //xx[SUM] -= 1;
        //if(xx[SUM] < 0) xx[SUM] += P;
        mk(x2,xx);
        conv(x1,x2);
        for(LL i = SUM; i < SUM + SUM + 1; i++)
            ans[i] -= x1[i];
        LL CNT = 0;
        for(LL i = SUM + 1; i < SUM + SUM + 1; i++)
            CNT += ans[i];
        ans[0] = 0;
        for(LL i = 1; i <= n; i ++)
            ans[0] += i *(n - i + 1);
        //cout << "CNT = " << CNT << " ans[0] = " << ans[0] << endl;
        ans[SUM] = ans[0] - CNT;
        for(LL i = SUM; i < SUM + SUM + 1; i++)
            printf("%I64d\n",ans[i]);
    }
}
```



  [1]: http://acm.hdu.edu.cn/showproblem.php?pid=1402
