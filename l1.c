/* ローレンツ方程式のカオスのストレンジアトラクタを描く.  */
/*
** ローレンツ方程式の説明についてはlorenz.hのコメント文を参照せよ. 
** アトラクタとは時間とともに体積が縮小する散逸力学系において
** 過渡状態経過後の定常状態のことである. 
** カオスのアトラクタを特に, ストレンジアトラクタト呼ぶ. 
** 常微分方程式の数値計算にはルンゲクッタ法を用いた. 
** ルンゲクッタ法の説明についてはrk.hのコメント文を参照せよ. 
** 描画にはgnuplotを用いた. 
** このプログラムはlinuxで動作確認している. 
*/


#include <stdio.h>
#include <stdlib.h>
#include "rk.h"
#include "lorenz.h"
 /* lorenz(double t,double X[],double dXdt[]) */
#define N 3 /* 状態変数の次元 */
#define h 0.01 /* ルンゲクッタステップ */
#define T 10000 /* ルンゲクッタの計算回数 */

main()
{
  int t,i;
  /* double X0[N],double X1[N]の配列を動的に作成 */
  double *X0=vector(N),*X1=vector(N);

  /* 初期値の設定 */
  X0[0]=10.0;
  X0[1]=20.0;
  X0[2]=30.0;

  /* メインパート */
  for(t=0;t<=T-1;t++)
    {
      for(i=0;i<=N-1;i++)
	{
	  printf("%f",X0[i]);
	  if(i==N-1)putchar('\n');
	  else putchar(' ');
	}
      /*
      ** 関数lorenzと時刻h*tのとき状態変数X0を渡し, 
      ** 時刻h*(t+1)のとき状態変数X1を得る. 
      */
      rk(h,N,lorenz,h*t,X0,X1);
      copy_vector(N,X1,X0);
    }

  return 0;
}
