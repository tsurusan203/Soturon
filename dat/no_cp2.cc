#include "noise_cancel.h"
#include "config.hh"
#include "fft.hh"
#include "mul.hh"

/*

  ノイズキャンセル処理
  
   入力
     xin: 1ブロック分のサブセンサ入力
     yin: 1ブロック分のメインセンサ入力
     
   出力
     eout: 処理結果

  
   1ブロック分のノイズキャンセル処理を行う
   過去FF個分の入力データを内部に保持する

*/
void noise_cancel(float eout[SHIFT], const short_v2 din[SHIFT])
{
#pragma HLS interface ap_fifo port=din    
#pragma HLS interface ap_fifo port=eout
  static float W_real[FF] = {0.0f};
  static float W_imag[FF] = {0.0f};
  static float ebuf[SHIFT] = {0.0f};
  float X_real[FF], X_imag[FF];
  float Ax[NPblk];
  float Px[NPblk];
  float D_real[NPblk], D_imag[NPblk];
  float Y_real[FF], Y_imag[FF];    
  float dtmp[FF];
  float etmp[S];
  float Ay[NPblk];
  float E_real[FF], E_imag[FF];
  float Ae[NPblk];
  float tmp[FF];
  float Px_ave[NPblk];
  float Ratio_ex[NPblk];
  float Ratio_yx[NPblk];
  float Etmp_real[NPblk], Etmp_imag[NPblk];
  float EE_real[NPblk], EE_imag[NPblk], EEE[NPblk];
  float dW_real[FF], dW_imag[FF];
  float dD_real[FF], dD_imag[FF];
  float e[SHIFT];
  float sum_E, sum_dD, gamma2;
  int i, j, k, ii;
  static float xbuff[FF] = {0.0f};
#pragma HLS BIND_STORAGE variable=xbuff type=ram_2p impl=bram  
  static float ytmp[FF] = {0.0f};
#pragma HLS BIND_STORAGE variable=ytmp type=ram_2p impl=bram    
  static int wp_x = 0, wp_y = 0;
  int rp_x, rp_y;


  /*
    
    事前計算しておくテーブル類
    
  */

  // 残差信号のオーバーラップ用と周波数分析用の窓
#include "WE_win.h"  
  
  // フィルタ修正時の低域と高域の抑制処理用
#include "Phi_win.h"  
  
  // パワー・相関計算用
#include "AVE_win.h"  
  
  // 非線形処理関連
#include "NLP.h"  
  
  // エコーの予測値の大きさ補正処理用
#include "wwin.h"  

  for (i = 0; i < SHIFT; i++) {
    // 入力信号を獲得
    const short_v2 d = din[i];
    const float x = d.x / (float) FRACTION_GAIN;
    const float y = d.y / (float) FRACTION_GAIN;
    
    // サブセンサ信号
    xbuff[wp_x] = x;
    wp_x = (wp_x + 1) % FF;
      
    // メインセンサ信号
    ytmp[wp_y] = y;
    wp_y = (wp_y + 1) % FF;
  }

  rp_x = wp_x;
  rp_y = (wp_y - (S + DELAY) + FF) % FF;

  // 全体の長さがFFとなるFFT・パワースペクトルの計算
  fft<float,float,float,float,float>(Ax, 1, Px, 1, X_real, X_imag, xbuff, 0, FF, rp_x, WE_win, 0);

  // サブセンサ信号とデジタルフィルタ係数の掛け算（周波数毎）して
  // 周波数領域のエコーの予測値を計算
  mul(D_real, D_imag, X_real, X_imag, W_real, W_imag);
  
  // Dを逆FFTして，時間信号へ
  ifft<float, float>(dtmp, D_real, D_imag);

  // 1ブロック分のエコキャン処理実行
  for (i = 0; i < S; i++) {
    etmp[i] = ytmp[(i + rp_y) % FF] - dtmp[L + i];
  }
  fft<float,float,float,float,float>(Ay, 1, Px, 0, Y_real, Y_imag, ytmp, L, S, rp_y, WE_win, 1);

  // 残差／サブセンサ信号のゲイン平均値  20200831 調和平均に変更
  gain_average_reg(Ratio_yx, Ax, Ay, AVE_win);    

  for (j = 0; j < LOOP; ++j) {
    // フィルタのリセット処理
    // 残差の前半にLサンプルの零データを詰め込んでからFFT
    // （残差信号のスペクトル）
    fft<float,float,float,float,float>(Ae, 1, Px, 0, E_real, E_imag, etmp, L, S, 0, WE_win, 1);      

    // MEMSセンサ信号より、残差信号が大きい周波数ビンのフィルタ係数は
    // ０にリセット
    for (i = 0; i < NPblk; i++) {
      if (Ay[i] < Ae[i] * 0.5f) {
        W_real[i] = 0.0f;
        W_imag[i] = 0.0f;
      }
    }

    // リセット処理した更新ベクトルを時間領域へ    
    ifft<float,float>(tmp, W_real, W_imag);

    // 時間領域で拘束をかけた更新ベクトルを再度周波数領域へ
    fft<float,float,float,float,float>(Ax, 0, Px, 0, W_real, W_imag, tmp, 0, L, 0, WE_win, 1);

    // リセット処理したフィルタでエコーの予測値を再計算
    mul(D_real, D_imag, X_real, X_imag, W_real, W_imag);

    // Dを逆FFTして，時間信号へ
    ifft<float,float>(dtmp, D_real, D_imag);
    
    // 1ブロック分のエコキャン処理の再実行
    for (i = 0; i < S; i++) {
      etmp[i] = ytmp[(i + rp_y) % FF] - dtmp[L + i];
    }

    // 残差の前半にLサンプルの零データを詰め込んでから
    // FFT<FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FFT_BUF_T>（残差信号のスペクトル）
    fft<float,float,float,float,float>(Ae, 1, Px, 0, E_real, E_imag, etmp, L, S, 0, WE_win, 1);

    // パワー・相関計算
    // サブセンサ信号のパワー平均値    
    conv_reg(Px_ave, Px, AVE_win);
    
    // 残差／サブセンサ信号のゲイン平均値  20200831 調和平均に変更
    gain_average_reg(Ratio_ex, Ax, Ae, AVE_win);
    
    for (i = 0; i < NPblk; i++) {
      float lambda, invPx, Phi;
      
      if (Ratio_ex[i] > Ratio_yx[i]) {
        Ratio_ex[i] = Ratio_yx[i];
      }      

      // 周波数毎の残差／サブセンサ信号比のリミッタ処理
      Phi = Ae[i] / (Ax[i] + EPSILON);
      if (Phi > Ratio_ex[i]) {
        Phi = Ratio_ex[i];
      }
      Phi = Phi * Phi_win[i];      

      // 更新式の分母に加える正則化パラメータλの動的計算
      lambda = Ae[i] * Ax[i] / (Phi + EPSILON) - Px[i];
      invPx = 1.0f / (Px[i] + lambda + Px_ave[i] + 0.1f);

      // フィルタ出力のレベルと実際の残差レベルがなるべく揃うように
      // フィルタ更新ベクトルを補正
      EE_real[i] = E_real[i] * invPx;
      EE_imag[i] = E_imag[i] * invPx;
      Etmp_real[i] = 0.0f;
      Etmp_imag[i] = 0.0f;      
      EEE[i] = Ae[i] * invPx * Px[i];
    }
    
    for (ii = 0; ii < FIDX2_num - 1; ii++) {
      float gain;
      float dot_Ae = 0.0f, dot_EEE = 0.0f;
      
      for (k = FIDX2[ii]; k < FIDX2[ii + 2] + 1; k++) {
        const float wEEE = EEE[k] * wwin[widx[ii] + k - FIDX2[ii]];
        dot_Ae  += Ae[k] * wEEE;
        dot_EEE += EEE[k] * wEEE;
      }
      // 周波数帯域毎にレベル補正ゲインを決定
      gain = dot_Ae / (dot_EEE + EPSILON);
      for (k = FIDX2[ii]; k < FIDX2[ii + 2] + 1; k++) {
        Etmp_real[k]
          = Etmp_real[k] + gain * EE_real[k] * wwin[widx[ii] + k - FIDX2[ii]];
        Etmp_imag[k]
          = Etmp_imag[k] + gain * EE_imag[k] * wwin[widx[ii] + k - FIDX2[ii]];
      }
    }
    
    for (i = 0; i < NPblk; i++) { 
      //周波数領域の更新ベクトル計算            
      dW_real[i] =   X_real[i] * Etmp_real[i] + X_imag[i] * Etmp_imag[i];
      dW_imag[i] = - X_imag[i] * Etmp_real[i] + X_real[i] * Etmp_imag[i]; 
    }    

    // 更新ベクトルを時間領域へ
    ifft<float,float>(tmp, dW_real, dW_imag);
    
    // 時間領域で拘束をかけた更新ベクトルを再度周波数領域へ        
    fft<float,float,float,float,float>(Ax, 0, Px, 0, dW_real, dW_imag, tmp, 0, L, 0, WE_win, 1);
    
    // 更新したフィルタで残差に残ったエコー成分を再度予測
    mul(dD_real, dD_imag, X_real, X_imag, dW_real, dW_imag);
    
    // 時間領域に戻して
    ifft<float, float>(dtmp, dD_real, dD_imag);
    
    // 残差に残ったエコー成分の予測値のスペクトル
    fft<float,float,float,float,float>(Ax, 0, Px, 0, dD_real, dD_imag, dtmp, L, S, L, WE_win, 1);

    // 残差に残ったエコー成分の予測値のスペクトルと実際の残差とのレベル合わせ
    // 重み付き周波数スペクトルの相関に変更
    sum_E  = 0.0f;
    sum_dD = 0.0f;
    for (i = 0; i < NPblk; i++) {
      float wdD_real = dD_real[i] * Phi_win[i];
      float wdD_imag = dD_imag[i] * Phi_win[i];
      sum_E  += ( E_real[i] * wdD_real +  E_imag[i] * wdD_imag);
      sum_dD += (dD_real[i] * wdD_real + dD_imag[i] * wdD_imag);
    }
    gamma2 = sum_E / (sum_dD + EPSILON);
    
    // フィルタ再更新
    for (i = 0; i < NPblk; i++) {
      W_real[i] = W_real[i] + mu * gamma2 * dW_real[i];
      W_imag[i] = W_imag[i] + mu * gamma2 * dW_imag[i];
    }

    // 再度エコキャン実行
    for (i = 0; i < S; i++) {
      etmp[i] = etmp[i] - mu * gamma2 * dtmp[L + i];
    }
  }
  
  // 非線形抑圧の実行
#ifdef NLP
  float avEh[FF], avEh2[FF], avEYh[FF], avEYh2[FF];
  float cmpE_real[NPblk], cmpE_imag[NPblk];

  fft<float,float,float,float,float>(Ae, 1, Px, 0, E_real, E_imag, etmp, L, S, 0, WE_win, 1);

  nlp_conv1_reg2(avEh,  Ae, AVE_win2);
  nlp_conv1_reg3(avEh2, Ae, AVE_win3);
  nlp_conv2_reg4(avEYh,  Ay, Ae, AVE_win4);
  nlp_conv2_reg5(avEYh2, Ay, Ae, AVE_win5);  
  
  for (i = 0; i < NPblk; i++) {
    const int k = i + (FF / 2 - 1);
    const float avEh_i  = NLP_winL[i] * avEh[k]  + NLP_winH[i] * avEh2[k];
    const float avEYh_i = NLP_winL[i] * avEYh[k] + NLP_winH[i] * avEYh2[k];
    float Eh, El;    
    
    if (Ae[i] > avEh_i) {
      Eh = Ae[i] - avEh_i;
      El = avEh_i;
    }
    else {
      Eh = 0.0f;
      El = Ae[i];
    }
    const float t = 1.0f - fminf(avEYh_i, EY_th) / EY_th;
    const float Att = t * t;
    const float cmp = 1.0f / (1.0f + NLP_depth * Att * Eh / (El + EPSILON));
    cmpE_real[i] = (El + cmp * Eh) * E_real[i] / (Ae[i] + EPSILON);
    cmpE_imag[i] = (El + cmp * Eh) * E_imag[i] / (Ae[i] + EPSILON);
  }

  ifft<float, float>(E_real, cmpE_real, cmpE_imag);

  for (i = 0; i < SHIFT; i++) {
    e[i] = ebuf[i] + E_real[L + i] * WE_win[L + i];
    ebuf[i] = E_real[L + SHIFT + i] * WE_win[L + SHIFT + i];
  }
#else
  for (i = 0; i < SHIFT; i++) {
    e[i] = ebuf[i] + etmp[i] * WE_win[L + i];
    ebuf[i] = etmp[SHIFT + i] * WE_win[L + SHIFT + i];
  }  
#endif
  for (i = 0; i < SHIFT; i++) {
      eout[i] = static_cast<float>(e[i]);
  }
}

