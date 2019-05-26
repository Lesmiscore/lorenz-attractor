/* ルンゲクッタ・ルーチン (rk.h) */
/*
** このヘッダファイルはN個の状態変数を持つ
** 1階の常微分方程式を解くためのものである. 
** 3次までのルンゲクッタ法を用いる. 
** コンピュータシュミレーションでは, 微分方程式
** dx(t)/dt=F(x(t),t)を解析的に求めることはできないので
** 時間tを離散化して求めることになる. 
** 4次のルンゲクッタ法を式で表すと以下のようになる. 
** d1=hF(x(t),t)
** d2=hF(x(t)+d1/2,t+h/2)
** d3=hF(x(t)+d2/2,t+h/2)
** d4=hF(x(t)+d3,t+h)
** x(t+h)=x(t)+d1/6+d2/3+d3/3+d4/6
** hは離散化の時間幅であるルンゲクッタステップである. 
*/

#include <stdlib.h>

double *vector(int N) /* ベクトル領域の確保 */
{
  /* サイズNの領域をmallocで動的に作成する */
  return (double *)malloc(N*sizeof(double));
}

void free_vector(double *v) /* ベクトル領域の解放 */
{
  /* ポインタvの領域を解放する */
  free(v);
}

void copy_vector(int N,double a[],double b[])
 /* ベクトルのコピー */
{
  int i;
  /* ベクトルaをベクトルbにコピーする */
  for(i=0;i<=N-1;i++)b[i]=a[i];
}

/*
** ルンゲクッタステップh,微分方程式の階数N,微分方程式が記述された
** 関数dXdt,時刻tのときの状態変数が入ったベクトルX0を渡し, 時刻t+
** h の時の状態変数X1を得る. dXdtはhoge(double t,double X[],doubl
** edXdt[])の型でユーザが定義する. 
*/
void rk(double h,int N,
	void (*dXdt)(double t,double X[],double dXdt[]),
	double t,double X0[],double X[]) /* ルンゲクッタ */
{
  int i;
  /* double d1[N],double d2[N],double d3[N]の配列を動的に作成 */
  double *d1=vector(N),*d2=vector(N),*d3=vector(N);
  /* double Xa[N]double,X[N]の配列を動的に作成 */
  double *Xa=vector(N),*dX=vector(N);

  /* d1=hF(x(t),t) */
  dXdt(t,X0,dX);
  for(i=0;i<=N-1;i++)
    {
      d1[i]=h*dX[i];
      Xa[i]=X0[i]+0.5*d1[i];
    }

  /* d2=hF(x(t)+d1/2,t+h/2) */
  dXdt(t+0.5*h,Xa,dX);
  for(i=0;i<=N-1;i++)
    {
      d2[i]=h*dX[i];
      Xa[i]=X0[i]+0.5*d2[i];
    }

  /* d3=hF(x(t)+d2/2,t+h/2) */ 
  dXdt(t+0.5*h,Xa,dX);
  for(i=0;i<=N-1;i++)
    {
      d3[i]=h*dX[i];
      Xa[i]=X0[i]+d3[i];
    }

  /* x(t+h)=x(t)+d1/6+d2/3+d3/3 */
  dXdt(t+h,Xa,dX);
  for(i=0;i<=N-1;i++)X[i]=X0[i]
  +(d1[i]+d2[i]*2+d3[i]*2+h*dX[i])/6.0;

  /* 各配列の領域の解放 */
  free_vector(d1);free_vector(d2);free_vector(d3);
  free_vector(Xa);free_vector(dX);
}
