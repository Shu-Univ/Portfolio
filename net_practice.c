#include "net_practice_head.h"  // 自作ヘッダは<>ではなく""で囲う

// 式の形は MATLAB のコードと同じ
// 数値は論文サプリメントの pdf に合わせていく

void make_cell_params(FILE* paramas_file, int n_neuron){
  
  if((paramas_file=fopen("net_practice_params.dat","w"))==NULL){
    printf("Can't open net_practice_params.dat File! \n");
  }else{
    for(int32_t i=0; i<n_neuron; i++){
      fprintf(paramas_file,"%f " ,0.55);
      fprintf(paramas_file,"%f " ,0.66);
      fprintf(paramas_file,"%f " ,1.0);
      fprintf(paramas_file,"%f\n",3.141592);
    }
  }
  fclose(paramas_file);
}

void initialize_cell_params(double **cell_params, FILE* params_file, int n_neuron){
  for(int32_t i=0; i<n_neuron; i++){
    fscanf(params_file,"%lf %lf %lf %lf\n",&cell_params[COORD_X][i],&cell_params[COORD_Y][i],&cell_params[DEND_SHAPE][i],&cell_params[DENDE_DIRECTION][i]);
  }
  fclose(params_file);
}
void output_env_params(int n_neuron){ // マクロ定数を書き出すためのファイル
  FILE* env_params=NULL;
  if((env_params=fopen("net_practice_env_params.dat","w"))==NULL){
    printf("Can't open Output or Input File! \n");
    exit(1);
  }
  fprintf(env_params,"%f %f %f\n",DEND_RANGE_CURL    ,DEGREE_CURL    ,JUNK_ABLE_AREA_CURL    );
  fprintf(env_params,"%f %f %f\n",DEND_RANGE_STRAIGHT,DEGREE_STRAIGHT,JUNK_ABLE_AREA_STRAIGHT);
  fprintf(env_params,"%d %d %d\n",N_RANDAM,n_neuron,100);
  fprintf(env_params,"%d %d %d"  ,CURL    ,STRAIGHT,100);
  fclose(env_params);
}

double dist_betwen_2_neu(double** cell_params, int i_1, int i_2){
  double dx = cell_params[COORD_X][i_1] - cell_params[COORD_X][i_2];
  double dy = cell_params[COORD_Y][i_1] - cell_params[COORD_Y][i_2];
  return sqrt(dx*dx + dy*dy);
}

double distance(double x_1, double y_1, double x_2, double y_2){
  double dx = x_1 - x_2;
  double dy = y_1 - y_2;
  return sqrt(dx*dx + dy*dy);
}

void make_randam_points(double randam_points[][4], double x_1, double y_1, double direct_1, double dend_shape_d){

  double rand_for_r, rand_for_theta;
  int dend_shape = (int)dend_shape_d;

  if(dend_shape==CURL){
    for(int32_t i=0; i<N_RANDAM; i++){
      rand_for_r     = drand48();
      rand_for_theta = drand48();
      randam_points[i][0] = x_1 + rand_for_r*DEND_RANGE_CURL*cos(direct_1 - DEGREE_CURL/2 + rand_for_theta*DEGREE_CURL);
      randam_points[i][1] = y_1 + rand_for_r*DEND_RANGE_CURL*sin(direct_1 - DEGREE_CURL/2 + rand_for_theta*DEGREE_CURL);
    }
  }else if(dend_shape==STRAIGHT){
    for (int32_t i=0; i<N_RANDAM/2; i++){
      rand_for_r     = drand48();
      rand_for_theta = drand48();
      randam_points[i][0] = x_1 + rand_for_r*DEND_RANGE_STRAIGHT*cos(direct_1 - DEGREE_STRAIGHT/2 + rand_for_theta*DEGREE_STRAIGHT);
      randam_points[i][1] = y_1 + rand_for_r*DEND_RANGE_STRAIGHT*sin(direct_1 - DEGREE_STRAIGHT/2 + rand_for_theta*DEGREE_STRAIGHT);
    }
    for (int32_t i=N_RANDAM/2; i<N_RANDAM; i++){
      rand_for_r     = drand48();
      rand_for_theta = drand48();
      randam_points[i][0] = x_1 + rand_for_r*DEND_RANGE_STRAIGHT*cos(direct_1+M_PI - DEGREE_STRAIGHT/2 + rand_for_theta*DEGREE_STRAIGHT);
      randam_points[i][1] = y_1 + rand_for_r*DEND_RANGE_STRAIGHT*sin(direct_1+M_PI - DEGREE_STRAIGHT/2 + rand_for_theta*DEGREE_STRAIGHT);
    }
  }
}

int count_included_points(double randam_points[][4],double x_2, double y_2, double direct_2, int start_index, int end_index, int dend_shape){
  
  int direction_type;
  int count = 0;
  double range  = 0.0;  // コンパイラがうるさいので仕方なく初期化
  double degree = 0.0;  // コンパイラがうるさいので仕方なく初期化
  double dist, rel_degree;
  double degree_start, degree_end;

  switch (dend_shape){
  case CURL:
    range  = DEND_RANGE_CURL;
    degree = DEGREE_CURL;
    break;
  case STRAIGHT:
    range  = DEND_RANGE_STRAIGHT;
    degree = DEGREE_STRAIGHT;
    break;
  default:
    printf("dend_shape = %d でエラーを起こしてます!!!（count_include_points 内）\n",dend_shape);
    break;
  }

  // 目標の dend が 0 を跨いでいるかどうかで場合分け
  degree_start = direct_2 - degree/2.0;
  degree_end   = direct_2 + degree/2.0;
  if((0.0 < degree_start)&&(degree_end < 2.0*M_PI)){
    direction_type = 1;
  }else if((degree_start < 0.0)||(2.0*M_PI < degree_end)){
    direction_type = -1;
    if(degree_start < 0.0){
      degree_start += 2.0*M_PI;
    }else if(2.0*M_PI < degree_end){
      degree_end   -= 2.0*M_PI;
    }
  }else{
    direction_type = 3;
  }
  printf("direction_type = %2d\n",direction_type);

  for(int32_t i=start_index; i<end_index; i++){
    dist = distance(randam_points[i][0],randam_points[i][1], x_2,y_2);
    if(dist < range){
      // 目標の cell に対する、ランダム座標の相対角度 rel_degree を作成
      rel_degree = atan2(randam_points[i][1]-y_2,randam_points[i][0]-x_2);
      if(rel_degree <= 0.0){
        rel_degree += 2.0 * M_PI; // 0~2π の範囲にする
      }

      // 目標 ION の樹状突起の角度内に入っていればカウント
      // 樹状突起が 0 を跨いでいるかどうかで判定条件が変わる
      if(direction_type == 1){
        if((degree_start < rel_degree)&&(rel_degree < degree_end)){
          count += 1;
        }
      }else if(direction_type == -1){
        if((rel_degree < degree_end)||(degree_start < rel_degree)){
          count += 1;
        }
      }else if(direction_type == 3){
        printf("エラーです！！！\n");
      }
    }
    randam_points[i][2] = dist;
    randam_points[i][3] = rel_degree;
  }
  return count;
}

void only_insert_connect_strength(double** connect, double** cell_params, double randam_points[N_RANDAM][4], int n_neuron){

  int dend_shape_2;
  int count_sum = 0;
  double overlapped_ratio, overlapped_area;

  for(int32_t i=0; i<n_neuron; i++){
    connect[i][i] = 0.0; //自分との結合は使われることはないが、0.0 を入れておく
  }

  for(int32_t i=0; i<n_neuron; i++){
    for(int32_t j=(i+1); j<(n_neuron); j++){
      
      // ここで、結合範囲内でランダム座標を生成するのは i 番目、それに対して目標セルになるのは j 番目とする
      make_randam_points(randam_points, cell_params[COORD_X][i], cell_params[COORD_Y][i], cell_params[DENDE_DIRECTION][i], cell_params[DEND_SHAPE][i]);
      dend_shape_2 = (int)cell_params[DEND_SHAPE][j];

      switch(dend_shape_2){
      case CURL:
        // count_included の最後の引数を cell_params[DEND_SHAPE][i] じゃなくて CURL にしてるけどいいのか？？
        // cell_params 内と dend_shape_2 と CURL は必ず同じになることは明白
        // その上 cell_params だと int へのキャストの一手間が入るので、int である CURL で渡した方がいい
        count_sum  = count_included_points(randam_points, cell_params[COORD_X][j],cell_params[COORD_Y][j],cell_params[DENDE_DIRECTION][j]       , 0, N_RANDAM, CURL);
        overlapped_ratio = (double)count_sum/(double)N_RANDAM;
        overlapped_area  = JUNK_ABLE_AREA_CURL     * overlapped_ratio;
        break;
      case STRAIGHT:
        count_sum  = count_included_points(randam_points, cell_params[COORD_X][j],cell_params[COORD_Y][j], cell_params[DENDE_DIRECTION][j]      , 0, N_RANDAM, STRAIGHT);
        count_sum += count_included_points(randam_points, cell_params[COORD_X][j],cell_params[COORD_Y][j],(cell_params[DENDE_DIRECTION][j]+M_PI), 0, N_RANDAM, STRAIGHT);
        overlapped_ratio = (double)count_sum/(double)N_RANDAM;
        overlapped_area  = JUNK_ABLE_AREA_STRAIGHT * overlapped_ratio;
        break;
      default:
        printf("dend_shape_2 = %d でエラーを起こしてます!!!（switch 内）\n",dend_shape_2);
        break;
      }
      printf("neuron_%d の結合可能面積の、neuron_%d に含まれる割合が %.2f %%なので、絶対面積は %.2f \n",i,j, overlapped_ratio*100.0, overlapped_area);
      connect[i][j] = overlapped_area * AREA_CONDACTANCE_RATIO;
      connect[j][i] = overlapped_area * AREA_CONDACTANCE_RATIO;
    }
  }
}

void calc_each_dvdt(double **vars, double** materials, double** connect, int n_neuron, int neuron_index, double t, double I_app){

  // V common
  double V_app     = 0.0;
  double dV_soma_dt, dV_dend_dt, dV_axon_dt;
  
  // ギャップ結合のため
  double v_deff, g_gap;
  double I_gap_each = 0.0;
  double I_gap = 0.0;

  // soma
  double dk_dt, dl_dt, dh_dt, dn_dt, dp_dt, dx_dt_soma;
  // dend
  double dq_dt, dr_dt, ds_dt, dCa_dt;
  // axon
  double dh_dt_axon, dx_dt_axon;
  

  //   0_soma   Update somatic   components
  dh_dt      = (inf_soma(H_INF     , vars[V_SOMA][neuron_index]) - materials[SODIUM_H][neuron_index]        )/tau_soma(TAU_H     , vars[V_SOMA][neuron_index]);
  dk_dt      = (inf_soma(K_INF     , vars[V_SOMA][neuron_index]) - materials[CALCIUM_K][neuron_index]       )/tau_soma(TAU_K     , vars[V_SOMA][neuron_index]);
  dl_dt      = (inf_soma(L_INF     , vars[V_SOMA][neuron_index]) - materials[CALCIUM_L][neuron_index]       )/tau_soma(TAU_L     , vars[V_SOMA][neuron_index]);
  dn_dt      = (inf_soma(N_INF     , vars[V_SOMA][neuron_index]) - materials[POTASSIUM_N][neuron_index]     )/tau_soma(TAU_N     , vars[V_SOMA][neuron_index]);
  dp_dt      = (inf_soma(P_INF     , vars[V_SOMA][neuron_index]) - materials[POTASSIUM_P][neuron_index]     )/tau_soma(TAU_N     , vars[V_SOMA][neuron_index]);
  dx_dt_soma = (inf_soma(X_INF_SOMA, vars[V_SOMA][neuron_index]) - materials[POTASSIUM_X_SOMA][neuron_index])/tau_soma(TAU_X_SOMA, vars[V_SOMA][neuron_index]);
  materials[SODIUM_H][neuron_index]         += DELTA * dh_dt          ;
  materials[CALCIUM_K][neuron_index]        += DELTA * dk_dt          ;
  materials[CALCIUM_L][neuron_index]        += DELTA * dl_dt          ;
  materials[M_SOMA][neuron_index]            = inf_soma(M_INF, vars[V_SOMA][neuron_index]);
  materials[POTASSIUM_N][neuron_index]      += DELTA * dn_dt          ;
  materials[POTASSIUM_P][neuron_index]      += DELTA * dp_dt          ;
  materials[POTASSIUM_X_SOMA][neuron_index] += DELTA * dx_dt_soma     ;

  //   0_dend   Update dendritic components
  dq_dt  = (inf_dend(Q_INF, vars[V_DEND][neuron_index])       - materials[HCURRENT_Q][neuron_index] )/tau_dend(TAU_Q, vars[V_DEND][neuron_index] );
  dr_dt  = (inf_dend(R_INF, vars[V_DEND][neuron_index])       - materials[CALCIUM_R][neuron_index]  )/tau_dend(TAU_R, vars[V_DEND][neuron_index] );
  ds_dt  = (inf_dend(S_INF, materials[CA2PLUS][neuron_index]) - materials[POTASSIUM_S][neuron_index])/tau_dend(TAU_S, materials[CA2PLUS][neuron_index]);
  dCa_dt = -3.0 * vars[CA_HIGH][neuron_index] - 0.075 * materials[CA2PLUS][neuron_index] * arbitrary;
  materials[HCURRENT_Q][neuron_index]  += DELTA * dq_dt ;
  materials[CALCIUM_R][neuron_index]   += DELTA * dr_dt ;
  materials[POTASSIUM_S][neuron_index] += DELTA * ds_dt ;
  materials[CA2PLUS][neuron_index]     += DELTA * dCa_dt;

  //   0_axon   Update axonal    components
  dh_dt_axon = (inf_axon(H_INF     , vars[V_AXON][neuron_index]) - materials[SODIUM_H_AXON][neuron_index]   )/tau_axon(TAU_H     , vars[V_AXON][neuron_index]);
  dx_dt_axon = (inf_axon(X_INF_AXON, vars[V_AXON][neuron_index]) - materials[POTASSIUM_X_AXON][neuron_index])/tau_axon(TAU_X_AXON, vars[V_AXON][neuron_index]);
  materials[SODIUM_H_AXON][neuron_index]    += DELTA * dh_dt_axon;
  materials[POTASSIUM_X_AXON][neuron_index] += DELTA * dx_dt_axon;
  materials[M_AXON][neuron_index]            = inf_axon(M_INF, vars[V_AXON][neuron_index]);

  //		1 Compute dendrite currents and update calcium
  // Soma-dendrite interaction current sd
  vars[SOMA_DEND][neuron_index] = (g_int / (1 - p1)) * (vars[V_DEND][neuron_index] - vars[V_SOMA][neuron_index]);
  // Inward high-threshold Ca current CaH
  vars[CA_HIGH  ][neuron_index] =  g_Ca_High   * materials[CALCIUM_R][neuron_index]   * materials[CALCIUM_R][neuron_index] * (vars[V_DEND][neuron_index] - V_Ca);
  // Outward Ca-dependent K current K_Ca
  vars[K_CA     ][neuron_index] =  g_K_Ca      * materials[POTASSIUM_S][neuron_index] * (vars[V_DEND][neuron_index] - V_K);
  // Leakage current ld
  vars[LEAK_DEND][neuron_index] =  g_leak_dend * (vars[V_DEND][neuron_index] - V_LEAK);
  // Inward anomalous rectifier h
  vars[H        ][neuron_index] =  g_h         * materials[HCURRENT_Q][neuron_index]  * (vars[V_DEND][neuron_index] - V_h);
  //		***** GABA A current *****
  vars[GAB_DEND ][neuron_index] = gbar_gaba_dend * g_gaba_dend * (vars[V_DEND][neuron_index] - V_gaba_dend);
  //		***** AMPA current *****
  vars[AMP_DEND ][neuron_index] = gbar_ampa_dend * g_ampa_dend * (vars[V_DEND][neuron_index] - V_ampa     );

  //		2 Compute somatic currents
  // Dendrite-soma interaction current
  vars[DEND_SOMA][neuron_index] = (g_int / p1) * (vars[V_SOMA][neuron_index] - vars[V_DEND][neuron_index] );
  // Inward low-threshold Ca current
  vars[CA_LOW   ][neuron_index] = g_Ca_Low * materials[CALCIUM_K][neuron_index] * materials[CALCIUM_K][neuron_index] * materials[CALCIUM_K][neuron_index] * materials[CALCIUM_L][neuron_index] * (vars[V_SOMA][neuron_index] - V_Ca);
  // Inward Na current
  vars[NA_SOMA  ][neuron_index] = g_Na_soma * materials[M_SOMA][neuron_index] * materials[M_SOMA][neuron_index] * materials[M_SOMA][neuron_index] * materials[SODIUM_H][neuron_index] * (vars[V_SOMA][neuron_index] - V_Na);
  // Leak current
  vars[LEAK_SOMA][neuron_index] = g_leak_soma * (vars[V_SOMA][neuron_index] - V_LEAK);
  // Potassium current
  vars[KDR_SOMA ][neuron_index] = g_Kdr_soma * materials[POTASSIUM_N][neuron_index] * materials[POTASSIUM_P][neuron_index]     * (vars[V_SOMA][neuron_index] - V_K);
  vars[K_SOMA   ][neuron_index] = g_K_soma   * pow(materials[POTASSIUM_X_SOMA][neuron_index], 4) * (vars[V_SOMA][neuron_index] - V_K);
  // Axon-soma interaction current
  vars[AXON_SOMA][neuron_index] = (g_int / (1 - p2)) * (vars[V_SOMA][neuron_index] - vars[V_AXON][neuron_index]);
  // ***** AMPA current *****
  vars[AMP      ][neuron_index] = gbar_ampa * g_ampa * (vars[V_DEND][neuron_index] - V_ampa); //	//	 これあってる ???????
  // ***** GABA A current *****
  vars[GAB_SOMA ][neuron_index] = gbar_gaba_soma * g_gaba_soma * (vars[V_SOMA][neuron_index] - V_gaba_soma);

  //		3 Compute axonic currents
  // Sodium
  vars[NA_AXON  ][neuron_index] = g_Na_axon    * materials[M_AXON][neuron_index] * materials[M_AXON][neuron_index] * materials[M_AXON][neuron_index] * materials[SODIUM_H_AXON][neuron_index] * (vars[V_AXON][neuron_index] - V_Na);
  // Leak
  vars[LEAK_AXON][neuron_index] = g_leak_axon  * (vars[V_AXON][neuron_index] - V_LEAK                    );
  // Soma-axon interaction current sa
  vars[SOMA_AXON][neuron_index] = (g_int / p2) * (vars[V_AXON][neuron_index] - vars[V_SOMA][neuron_index]);
  // Potassium (transient)
  vars[K_AXON   ][neuron_index] = g_K_axon     * pow(materials[POTASSIUM_X_AXON][neuron_index], 4) * (vars[V_AXON][neuron_index] - V_K);
  
  //   4 Compute leaks from soma-dendrite and soma-axon 
  // update voltages
  //dV_soma_dt = (-(I_dend_soma + I_axon_soma + I_leak_soma ))/Cm;
  //dV_dend_dt = (-(I_soma_dend               + I_leak_dend + I_cx36) + I_app)/Cm;
  //dV_axon_dt = (-(I_K_axon + I_soma_axon               + I_leak_axon + I_Na_axon))/Cm;

  // neuron_index が 1 の時のみ、SPIKE_TIMIMG の周辺時間のみ入力入れてみる
  if((neuron_index != 1)||(t/DELTA < SPIKE_TIMING*100)||((SPIKE_TIMING+DURATION)*100 < t/DELTA)){
    I_app = 0.0;
  }

  // 他の全ての樹状突起について、今の ION とのギャップ電流を加算
  for(int32_t i=0; i<n_neuron; i++){
    v_deff = (vars[V_DEND][i] - vars[V_DEND][neuron_index]);
    //g_gap = 0.8 * exp(-(v_deff * v_deff)/100.0) + 0.2;
    g_gap = connect[neuron_index][i] * (0.8 * exp(-(v_deff * v_deff)/100.0) + 0.2);
    I_gap_each = g_gap * v_deff;   // コンダクタンス × 電位差（V_[i] - V[index]）
    I_gap += I_gap_each;
  }
  
  dV_soma_dt = (-(vars[CA_LOW][neuron_index]  + vars[DEND_SOMA][neuron_index] + vars[AXON_SOMA][neuron_index] + vars[LEAK_SOMA][neuron_index] + vars[NA_SOMA][neuron_index] + vars[KDR_SOMA][neuron_index] + vars[K_SOMA][neuron_index] ))/Cm;
  dV_dend_dt = (-(vars[CA_HIGH][neuron_index] + vars[SOMA_DEND][neuron_index]                                 + vars[LEAK_DEND][neuron_index] + vars[K_CA][neuron_index]    + vars[H][neuron_index]) + I_app + I_gap)/Cm;
  dV_axon_dt = (-(vars[K_AXON][neuron_index]  + vars[SOMA_AXON][neuron_index]                                 + vars[LEAK_AXON][neuron_index] + vars[NA_AXON][neuron_index]))/Cm;
    
  //   5 Integrate dV/dt
  if(V_app == 0){
    vars[V_SOMA][neuron_index] += DELTA * dV_soma_dt;
  }else{
    vars[V_SOMA][neuron_index] = V_app;
  }
  vars[V_DEND][neuron_index] += DELTA * dV_dend_dt;
  vars[V_AXON][neuron_index] += DELTA * dV_axon_dt;

}



int main(int argc, char *argv[]){
  // argc はCL引数の個数を表す。argv 配列の中にコマンドライン引数の値たちが入る。
  // ただし、実行するプログラム名自体もCL引数に含まれるので、argv[0] = 実行ファイル名 になる

  print_hello(1);
  
  int n_neuron     = atoi(argv[1]);
  int neuron_index = atoi(argv[2]);
  double I_app     = atof(argv[3]);
  
  // 電流と電圧の各値の要素
  double** vars              = make_malloc_matrix(N_VARS     , n_neuron);
  // コンダクタンスを算出するためのゲート変数などの要素
  double** materials         = make_malloc_matrix(N_MATERIALS, n_neuron);
  // ION の各パラメータを格納する配列
  double** cell_params       = make_malloc_matrix(N_PARAMS   , n_neuron);
  // ION 同士の結合の有無を示す配列
  double** connect_strenghth = make_malloc_matrix(n_neuron   , n_neuron);
  // モンテカルロ用のランダム座標
  double randam_points[N_RANDAM][4];

  // 試しに初期化してみる
  for (int32_t i=0; i<n_neuron; i++){
    for (int32_t j=0; j<n_neuron; j++){
      connect_strenghth[i][j] = 0.0;
    }
  }

  FILE* params_file=NULL;
  FILE* fo=NULL;
  int output_place, output_style;
  output_place = DAT_OUTPUT;
  if(neuron_index<=-1){
    output_style = MULTI_CELL_ONLY_V;
  }else{ 
    output_style = SINGLE_CELL_I_V;
  }
  // ファイル開けるかの確認は、display 関数の中でやってまうと、開くたびに中身がリセットされるので、ここで一度だけ開く
  if(output_place==DAT_OUTPUT){
    if(((params_file=fopen("net_practice_params.dat","r"))==NULL)||((fo=fopen("net_practice.dat","w"))==NULL)){
      printf("Can't open Output or Input File! \n");
      exit(1);
    }
  }
  // モンテカルロ用
  time_t tp;
  time(&tp);
  srand48(tp);

  //make_cell_params(f_params,n_neuron);

  // initialize_vars の内部で for しないのは、好きな時にセル指定で初期化できる関数として使うため
  for(int32_t i=0; i<n_neuron; i++){
    initialize_vars(vars, materials, i);
  }

  initialize_cell_params(cell_params, params_file, n_neuron);

  vars_display(vars, neuron_index, n_neuron, 0.0, output_place, output_style, fo);

  only_insert_connect_strength(connect_strenghth, cell_params, randam_points, n_neuron);

  //単位時間の数だけループ
  for(int32_t nt = 1; nt < NT; nt++){
    double t = DELTA * nt;
    
    // 細胞の数だけループ
    for(int32_t i=0; i<n_neuron; i++){
      calc_each_dvdt(vars, materials, connect_strenghth, n_neuron, i, t, I_app);
    }

    if(t<70){// index=1 だけスタートを 70 だけ遅らせる
      initialize_vars(vars, materials,1);
    }

    vars_display(vars, neuron_index, n_neuron, t, output_place, output_style, fo);
  }

  output_env_params(n_neuron); // Python で描画する用

  printf("\ncell_params = \n");
  for(int32_t i=0; i<4; i++){
    for(int32_t j=0; j<n_neuron; j++){
      printf(" %5.1f ",cell_params[i][j]);
    }
    printf(" \n");
  }
  printf(" \n");
  printf("connect_strength = \n\n");
  for(int32_t i=0; i<n_neuron; i++){
    printf(" ");
    for(int32_t j=0; j<n_neuron; j++){
      printf("%5.2f ",connect_strenghth[i][j]);
      if(j%5==0){printf("  ");}
    }
    printf(" \n");
    if(i%5==0){printf("\n");}
  }

  if(output_style==DAT_OUTPUT){
    fclose(fo);
  }

  return 0;
}