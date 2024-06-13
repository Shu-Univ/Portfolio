#include "net_practice_head.h"

void print_hello(int n){
  for(int32_t i=0; i<n; i++){
  // make ファイルを更新したときは、反映されてるか確認のために、ここの文章を少し変えること
  printf("自作ヘッダの正常動作の確認! \n");
  }
}

//alpha や beta の計算は、モノによって形が大きく違うのでそれぞれコードする
double alpha_x_soma(const double V_soma){return (0.13 * (V_soma + 25.0)) / (1.0 - exp( -(V_soma + 25.0)/10.0 ));}
double  beta_x_soma(const double V_soma){return 1.69 * exp( -0.0125*(V_soma + 35.0) );}
double alpha_x_axon(const double V_axon){return (0.13 * (V_axon + 25.0)) / (1.0 - exp( -(V_axon + 25.0)/10.0 ));}
double  beta_x_axon(const double V_axon){return 1.69 * exp( -0.0125*(V_axon + 35.0) );}
//double alpha_soma(void){return min([0.00002 * Ca2Plus, 0.01]);}
double alpha_soma(const double Ca2Plus){return (0.00002*Ca2Plus >= 0.01)?0.01:0.00002*Ca2Plus;}
double  beta_soma(void){return 0.015;}
double alpha_r(const double V_dend){return 1.7 / (1.0 + exp( (-(V_dend - 5.0) )/ 13.9));}
double  beta_r(const double V_dend){return 0.02 * (V_dend + 8.5) / (exp((V_dend + 8.5)/5.0) - 1.0);}

double tau_soma(const int var, const double V_soma){      
  switch(var){
  case TAU_K:
    return 1.0;
  case TAU_L:
    return    ((20.0 * exp( (V_soma+160.0)/30.0 ))/(1.0+exp( (V_soma+84.0)/7.3 )))+35.0;
  case TAU_H:
    return         3.0 * exp( (-40.0 - V_soma)/33.0  );
  case TAU_N:
    return 5.0 + (47.0 * exp(-(-50.0 - V_soma)/900.0));
  case TAU_P:
    return ( 1.0 /(1.0 + exp((    -V_soma   - 51.0)/ -12.0))  );   // TAU_N と同じでいいんじゃないの？
  case TAU_X_SOMA:
    return 1.0/( alpha_x_soma(V_soma) + beta_x_soma(V_soma) );
  default:
    fprintf (stderr, "プログラム内 Switch で Error: no such vars: %d\n", var);
    exit(1);
  }
}
double tau_dend(const int var, const double V_dend){
  switch(var){
  case TAU_Q:
    return 1.0 /(exp(-0.086 * V_dend - 14.6) + exp(0.070 * V_dend - 1.87));
  case TAU_S:
    // ここだけ、V_dend という名の Ca2Plus の値を入れる
    return 1.0 / (alpha_soma(V_dend) + beta_soma()   );
  case TAU_R:
    return 5.0 / (alpha_r(V_dend)    + beta_r(V_dend));
  default:
    fprintf (stderr, "プログラム内 Switch で Error: no such vars: %d\n", var);
    exit(1);
  }
}
double tau_axon(const int var, const double V_axon){
  switch(var){
  case TAU_H:
    return 1.5 * exp( (-40.0 - V_axon)/33.0 );
  case TAU_X_AXON:
    return 1.0/( alpha_x_axon(V_axon) + beta_x_axon(V_axon) );
  default:
    fprintf (stderr, "プログラム内 Switch文 で Error: no such vars: %d\n", var);
    exit(1);
  }
}
double inf_soma(const int var, const double V_soma){
  switch(var){
  case K_INF:
    return (1.0/( 1.0 + exp( -1.0 * (V_soma + 61.0 )/4.2 ) ));
  case L_INF:
    return (1.0/( 1.0 + exp( (    V_soma + 85.5)/8.5 )));
  case M_INF:
    return  1.0/( 1.0 + exp( (-30.0 -V_soma       )/ 5.5 ));
  case H_INF:
    return  1.0/( 1.0 + exp( (-70.0 -V_soma       )/-5.8 ));
  case N_INF:
    return  1.0/( 1.0 + exp( (-3.0  -V_soma       )/10.0 ));
  case P_INF:
    return  1.0/( 1.0 + exp( (-V_soma       -51.0)/-12.0 ) );
  case X_INF_SOMA:
    return alpha_x_soma(V_soma)/( alpha_x_soma(V_soma) + beta_x_soma(V_soma) );
  default:
    fprintf (stderr, "プログラム内 Switch で Error: no such vars: %d\n", var);
    exit(1);
  }
}
double inf_dend(const int var, const double V_dend){
  switch(var){
  case Q_INF:
    return 1.0/( 1.0 + exp((V_dend + 80.0) / 4.0) );
  case R_INF:
    return alpha_r(V_dend)/( alpha_r(V_dend) + beta_r(V_dend) );
  case S_INF:
    // ここだけ、V_dend という名の Ca2Plus の値を入れる
    return alpha_soma(V_dend)/(alpha_soma(V_dend) + beta_soma());
  default:
    fprintf (stderr, "プログラム内 Switch で Error: no such vars: %d\n", var);
    exit(1);
  }
}
double inf_axon(const int var, const double V_axon){
  switch(var){
  case M_INF:
    return 1.0/(1.0 + exp( (-30.0 - V_axon)/ 5.5   ) );
  case H_INF:
    return 1.0/(1.0 + exp( (-60.0 - V_axon)/(-5.8) ) );
  case X_INF_AXON:
    return alpha_x_axon(V_axon)/( alpha_x_axon(V_axon) + beta_x_axon(V_axon) );
  default:
    fprintf (stderr, "プログラム内 Switch で Error: no such vars: %d\n", var);
    exit(1);
  }
}


double** make_malloc_matrix(int low_size, int col_size){
  double** matrix_addres = (double**)malloc(sizeof(double*) * low_size);
  for(int32_t i=0; i<low_size; i++){
    matrix_addres[i] = (double*)malloc(sizeof(double*) * col_size);
  }
  return matrix_addres;
}


void initialize_vars(double** vars, double** materials, int neuron_index){

  // soma materials
  materials[SODIUM_H        ][neuron_index] = inf_soma(H_INF     , V_SOMA_INIT);
  materials[CALCIUM_K       ][neuron_index] = inf_soma(K_INF     , V_SOMA_INIT);
  materials[CALCIUM_L       ][neuron_index] = inf_soma(L_INF     , V_SOMA_INIT);
  materials[M_SOMA          ][neuron_index] = inf_soma(M_INF     , V_SOMA_INIT);
  materials[POTASSIUM_N     ][neuron_index] = inf_soma(N_INF     , V_SOMA_INIT);
  materials[POTASSIUM_P     ][neuron_index] = inf_soma(P_INF     , V_SOMA_INIT);
  materials[POTASSIUM_X_SOMA][neuron_index] = inf_soma(X_INF_SOMA, V_SOMA_INIT);
  // dend materials
  materials[HCURRENT_Q      ][neuron_index] = inf_dend(Q_INF     , V_DEND_INIT);
  materials[CALCIUM_R       ][neuron_index] = inf_dend(R_INF     , V_DEND_INIT);
  // axon materials
  materials[SODIUM_H_AXON   ][neuron_index] = inf_axon(H_INF     , V_AXON_INIT);
  materials[POTASSIUM_X_AXON][neuron_index] = inf_axon(X_INF_AXON, V_AXON_INIT);
  materials[M_AXON          ][neuron_index] = inf_axon(M_INF     , V_AXON_INIT);

  // ここから電位の初期化
  vars[V_DEND][neuron_index] = V_DEND_INIT;
  vars[V_SOMA][neuron_index] = V_SOMA_INIT;
  vars[V_AXON][neuron_index] = V_AXON_INIT;

  // ここから電流の初期化
  // ここから dend
  vars[SOMA_DEND][neuron_index] = (g_int / (1 - p1)) * 0.0;
  vars[CA_HIGH  ][neuron_index] = g_Ca_High * materials[CALCIUM_R][neuron_index] * materials[CALCIUM_R][neuron_index] * (V_DEND_INIT - V_Ca);
  vars[LEAK_DEND][neuron_index] = g_leak_dend * (V_DEND_INIT - V_LEAK);
  vars[H        ][neuron_index] = g_h * materials[HCURRENT_Q][neuron_index] * (V_DEND_INIT - V_h);
  vars[GAB_DEND ][neuron_index] = 0.0;
  vars[AMP_DEND ][neuron_index] = 0.0;

  // ここから soma
  vars[DEND_SOMA][neuron_index] = (g_int / p1)* (V_SOMA_INIT - V_DEND_INIT );
  vars[CA_LOW   ][neuron_index] = g_Ca_Low    * materials[CALCIUM_K][neuron_index] * materials[CALCIUM_K][neuron_index] * materials[CALCIUM_K][neuron_index] * materials[CALCIUM_L][neuron_index] * (V_SOMA_INIT - V_Ca);
  vars[NA_SOMA  ][neuron_index] = g_Na_soma   * materials[M_SOMA][neuron_index] * materials[M_SOMA][neuron_index] * materials[M_SOMA][neuron_index] * materials[SODIUM_H][neuron_index] * (V_SOMA_INIT - V_Na);
  vars[LEAK_SOMA][neuron_index] = g_leak_soma * (V_SOMA_INIT - V_LEAK);
  vars[KDR_SOMA ][neuron_index] = g_Kdr_soma  * materials[POTASSIUM_N][neuron_index] * materials[POTASSIUM_P][neuron_index] * (V_SOMA_INIT - V_K);
  vars[K_SOMA   ][neuron_index] = g_K_soma    * materials[POTASSIUM_X_SOMA][neuron_index] * materials[POTASSIUM_X_SOMA][neuron_index] * materials[POTASSIUM_X_SOMA][neuron_index] * materials[POTASSIUM_X_SOMA][neuron_index] * (V_SOMA_INIT - V_K);
  vars[AXON_SOMA][neuron_index] = (g_int / (1 - p2)) * (V_SOMA_INIT - V_AXON_INIT);
  vars[AMP      ][neuron_index] = 0.0;
  vars[GAB_SOMA ][neuron_index] = 0.0;

  // ここから axon
  vars[NA_AXON  ][neuron_index] = g_Na_axon    * materials[M_AXON][neuron_index] * materials[M_AXON][neuron_index] * materials[M_AXON][neuron_index] * materials[SODIUM_H_AXON][neuron_index] * (V_AXON_INIT - V_Na);
  vars[LEAK_AXON][neuron_index] = g_leak_axon  * (V_AXON_INIT - V_LEAK);
  vars[SOMA_AXON][neuron_index] = (g_int / p2) * (V_AXON_INIT - V_SOMA_INIT);
  vars[K_AXON   ][neuron_index] = g_K_axon * materials[POTASSIUM_X_AXON][neuron_index] * materials[POTASSIUM_X_AXON][neuron_index] * materials[POTASSIUM_X_AXON][neuron_index] * materials[POTASSIUM_X_AXON][neuron_index] * (V_AXON_INIT - V_K);
  // この３つの初期化は、電流の初期化が終わってからする
  materials[CA2PLUS    ][neuron_index] = -1.0 * vars[CA_HIGH][neuron_index] * 3.0 / 0.075	;
  materials[POTASSIUM_S][neuron_index] = inf_dend(S_INF, materials[CA2PLUS][neuron_index]);
  vars[K_CA][neuron_index]  = g_K_Ca * materials[POTASSIUM_S][neuron_index] * (vars[V_DEND][neuron_index] - V_K);
}

void vars_display(double** vars, int neuron_index, int neuron_n, double t, int output_place, int output_style, FILE* fo){

  switch(output_place){
  case STD_OUTPUT:
    printf("%f "  ,t);
    switch(output_style){
    case SINGLE_CELL_I_V:
      printf("%f "  ,t);
      printf("%f "  ,vars[V_SOMA][neuron_index]);
      printf("%f "  ,vars[V_DEND][neuron_index]);
      printf("%f "  ,vars[V_AXON][neuron_index]);

      //printf("%f "  ,vars[I_SOMA_DEND][neuron_index]);  
      //printf("%f "  ,vars[I_CA_HIGH  ][neuron_index]);
      //printf("%f "  ,vars[I_K_CA     ][neuron_index]);
      //printf("%f "  ,vars[I_LEAK_DEND][neuron_index]);
      //printf("%f "  ,vars[I_H        ][neuron_index]);
      //printf("%f "  ,vars[I_GAB_DEND ][neuron_index]);
      //printf("%f "  ,vars[I_AMP_DEND ][neuron_index]);

      //printf("%f "  ,vars[I_DEND_SOMA][neuron_index]);
      //printf("%f "  ,vars[I_CA_LOW   ][neuron_index]);
      //printf("%f "  ,vars[I_NA_SOMA  ][neuron_index]);
      //printf("%f "  ,vars[I_LEAK_SOMA][neuron_index]);
      //printf("%f "  ,vars[I_KDR_SOMA ][neuron_index]);
      //printf("%f "  ,vars[I_K_SOMA   ][neuron_index]);
      //printf("%f "  ,vars[I_AXON_SOMA][neuron_index]);
      //printf("%f "  ,vars[I_AMP      ][neuron_index]);
      //printf("%f "  ,vars[I_GAB_SOMA ][neuron_index]);

      //printf("%f "  ,vars[I_NA_AXON  ][neuron_index]);
      //printf("%f "  ,vars[I_LEAK_AXON][neuron_index]);
      //printf("%f "  ,vars[I_SOMA_AXON][neuron_index]);
      //printf("%f " 	,vars[I_K_AXON   ][neuron_index]);
      break;
    case MULTI_CELL_ONLY_V:
      for (int32_t i=0; i<neuron_n; i++){
        printf("%f "  ,vars[V_SOMA][i]);
        printf("%f "  ,vars[V_DEND][i]);
        printf("%f "  ,vars[V_AXON][i]);
      }
      break;
    }
    printf("\n");
    break;

  case DAT_OUTPUT:
    fprintf(fo, "%f " ,t);

    switch(output_style){
    case SINGLE_CELL_I_V:
    fprintf(fo, "%f " ,vars[V_SOMA][neuron_index]);
    fprintf(fo, "%f " ,vars[V_DEND][neuron_index]);
    fprintf(fo, "%f " ,vars[V_AXON][neuron_index]);
    fprintf(fo, "%f " ,vars[SOMA_DEND][neuron_index]);  
    fprintf(fo, "%f " ,vars[CA_HIGH  ][neuron_index]);
    fprintf(fo, "%f " ,vars[K_CA     ][neuron_index]);
    fprintf(fo, "%f " ,vars[LEAK_DEND][neuron_index]);
    fprintf(fo, "%f " ,vars[H        ][neuron_index]);
    fprintf(fo, "%f " ,vars[GAB_DEND ][neuron_index]);
    fprintf(fo, "%f " ,vars[AMP_DEND ][neuron_index]);

    fprintf(fo, "%f " ,vars[DEND_SOMA][neuron_index]);
    fprintf(fo, "%f " ,vars[CA_LOW   ][neuron_index]);
    fprintf(fo, "%f " ,vars[NA_SOMA  ][neuron_index]);
    fprintf(fo, "%f " ,vars[LEAK_SOMA][neuron_index]);
    fprintf(fo, "%f " ,vars[KDR_SOMA ][neuron_index]);
    fprintf(fo, "%f " ,vars[K_SOMA   ][neuron_index]);
    fprintf(fo, "%f " ,vars[AXON_SOMA][neuron_index]);
    fprintf(fo, "%f " ,vars[AMP      ][neuron_index]);
    fprintf(fo, "%f " ,vars[GAB_SOMA ][neuron_index]);

    fprintf(fo, "%f "  ,vars[NA_AXON  ][neuron_index]);
    fprintf(fo, "%f "  ,vars[LEAK_AXON][neuron_index]);
    fprintf(fo, "%f "  ,vars[SOMA_AXON][neuron_index]);
    fprintf(fo, "%f "  ,vars[K_AXON   ][neuron_index]);
      break;
    case MULTI_CELL_ONLY_V:
      for(int32_t i=0; i<neuron_n; i++){
      fprintf(fo, "%f " ,vars[V_SOMA][i]);
      fprintf(fo, "%f " ,vars[V_DEND][i]);
      fprintf(fo, "%f " ,vars[V_AXON][i]);
      }
      break;
    }
    fprintf(fo,"\n");
    break;
  }
}

