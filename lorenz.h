/* ローレンツ方程式 (lorenz.h) */
/*
** このヘッダファイルはローレンツ方程式の関数の型を
** ルンゲクッタ・ルーチン (rk.h)のために記述する. 
** ローレンツ方程式とは, 大気の対流現象の大幅な近似の下でのモデルで, 
** 3変数の常微分方程式である. 
** dx/dt=-σx+σy
** dy/dt=-xz+rx-y
** dz/dt=xy-bz
** σ,b,rはパラメータであり, ここでは
** σ=10,b=8/3,r=28とする. 
** ローレンツ方程式は初期値に対して鋭敏に反応する
** これは, バタフライ効果の語源になった現象である. 
*/

/* 関数および変数の表記をrk.hと整合させる */
#define dxdt dXdt[0]
#define dydt dXdt[1]
#define dzdt dXdt[2]
#define x X[0]
#define y X[1]
#define z X[2]

/* ローレンツ方程式 */
void lorenz(double t,double X[],double dXdt[])
{
  double sigma=10.0,b=8.0/3,r=28.0;
  dxdt=-sigma*x+sigma*y;
  dydt=-x*z+r*x-y;
  dzdt=x*y-b*z;
}

/* 定義の解除 */
#undef dxdt
#undef dydt
#undef dzdt
#undef x
#undef y
#undef z
