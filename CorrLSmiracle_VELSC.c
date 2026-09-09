/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CorrLSmiracle_VELSC.c
 *
 * Code generation for function 'CorrLSmiracle_VELSC'
 *
 */

/* Include files */
#include "CorrLSmiracle_VELSC.h"
#include "CorrLSmiracle_VELSC_emxutil.h"
#include "CorrLSmiracle_VELSC_types.h"
#include "correfMcss.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Declarations */
static void binary_expand_op(emxArray_real_T *in1, const emxArray_real_T *in2,
                             const emxArray_real_T *in3);

/* Function Definitions */
static void binary_expand_op(emxArray_real_T *in1, const emxArray_real_T *in2,
                             const emxArray_real_T *in3)
{
  emxArray_real_T *b_in1;
  const double *in2_data;
  const double *in3_data;
  double *b_in1_data;
  double *in1_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  int stride_2_1;
  in3_data = in3->data;
  in2_data = in2->data;
  in1_data = in1->data;
  emxInit_real_T(&b_in1, 2);
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = 1;
  if (in3->size[1] == 1) {
    if (in2->size[1] == 1) {
      loop_ub = in1->size[1];
    } else {
      loop_ub = in2->size[1];
    }
  } else {
    loop_ub = in3->size[1];
  }
  b_in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(b_in1, i);
  b_in1_data = b_in1->data;
  stride_0_1 = (in1->size[1] != 1);
  stride_1_1 = (in2->size[1] != 1);
  stride_2_1 = (in3->size[1] != 1);
  for (i = 0; i < loop_ub; i++) {
    b_in1_data[i] = (in1_data[i * stride_0_1] + in2_data[i * stride_1_1]) +
                    in3_data[i * stride_2_1];
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(in1, i);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = b_in1_data[i];
  }
  emxFree_real_T(&b_in1);
}

void CorrLSmiracle_VELSC(const double g[1151], double seed, double N, double Lo,
                         double K, double Code, emxArray_real_T *Ry_gk)
{
  static const signed char H[16] = {1, 1, 1,  1,  1, -1, 1,  -1,
                                    1, 1, -1, -1, 1, -1, -1, 1};
  static const signed char iv[4] = {1, 1, 2, 2};
  emxArray_real_T *C;
  emxArray_real_T *S;
  emxArray_real_T *b_y0;
  emxArray_real_T *setL0;
  emxArray_real_T *setL1;
  emxArray_real_T *x0;
  double d;
  double *C_data;
  double *S_data;
  double *setL0_data;
  double *setL1_data;
  double *x0_data;
  int b_i;
  int b_loop_ub_tmp;
  int i;
  int i1;
  int i2;
  int loop_ub_tmp;
  /*  CorrLS_modificado: Implementa el correlador adaptado al segundo diagrama
   */
  /*  Entradas: */
  /*    g: Señal de entrada */
  /*    seed: Semilla para secuencia complementaria */
  /*    N: Número de secuencias CCS */
  /*    Lo: Longitud de secuencias CCS */
  /*    K: Número de códigos */
  /*    Code: Índice de fila en matriz de Hadamard */
  /*  Parámetros iniciales */
  /*  Golay */
  /*  Longitud de zona libre de interferencia */
  /*  Mitad de los códigos */
  /*  Permutación */
  /* Matriz de Hadamard  */
  /*  Obtener códigos complementarios desde las semillas */
  /*  BLOQUE 1 2: correlación con códigos complementarios */
  emxInit_real_T(&x0, 2);
  i = x0->size[0] * x0->size[1];
  x0->size[0] = 1;
  x0->size[1] = (int)Lo + 1151;
  emxEnsureCapacity_real_T(x0, i);
  x0_data = x0->data;
  for (i = 0; i < 1151; i++) {
    x0_data[i] = g[i];
  }
  loop_ub_tmp = (int)Lo;
  for (i = 0; i < loop_ub_tmp; i++) {
    x0_data[i + 1151] = 0.0;
  }
  emxInit_real_T(&setL0, 2);
  correfMcss(N, x0, Lo, seed, setL0);
  setL0_data = setL0->data;
  i = x0->size[0] * x0->size[1];
  x0->size[0] = 1;
  x0->size[1] = (int)Lo + 1151;
  emxEnsureCapacity_real_T(x0, i);
  x0_data = x0->data;
  for (i = 0; i < 1151; i++) {
    x0_data[i] = g[i];
  }
  for (i = 0; i < loop_ub_tmp; i++) {
    x0_data[i + 1151] = 0.0;
  }
  emxInit_real_T(&setL1, 2);
  if (seed < Lo / 2.0) {
    d = Lo / 2.0 + seed;
  } else {
    d = seed - Lo / 2.0;
  }
  correfMcss(N, x0, Lo, d, setL1);
  setL1_data = setL1->data;
  /*  Separación de componentes C y S */
  /*  <-- CAMBIO MÍNIMO (antes era 2*Lo) */
  emxInit_real_T(&C, 2);
  i = C->size[0] * C->size[1];
  C->size[0] = 2;
  loop_ub_tmp = (int)(Lo + 1151.0);
  C->size[1] = (int)(Lo + 1151.0);
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  b_loop_ub_tmp = (int)(Lo + 1151.0) << 1;
  for (i = 0; i < b_loop_ub_tmp; i++) {
    C_data[i] = 0.0;
  }
  emxInit_real_T(&S, 2);
  i = S->size[0] * S->size[1];
  S->size[0] = 2;
  S->size[1] = (int)(Lo + 1151.0);
  emxEnsureCapacity_real_T(S, i);
  S_data = S->data;
  for (i = 0; i < b_loop_ub_tmp; i++) {
    S_data[i] = 0.0;
  }
  for (i = 0; i < loop_ub_tmp; i++) {
    C_data[2 * i] = setL0_data[setL0->size[0] * i];
  }
  for (i = 0; i < loop_ub_tmp; i++) {
    C_data[2 * i + 1] = setL1_data[setL1->size[0] * i];
  }
  for (i = 0; i < loop_ub_tmp; i++) {
    S_data[2 * i] = setL0_data[setL0->size[0] * i + 1];
  }
  emxFree_real_T(&setL0);
  for (i = 0; i < loop_ub_tmp; i++) {
    S_data[2 * i + 1] = setL1_data[setL1->size[0] * i + 1];
  }
  emxFree_real_T(&setL1);
  /*  Lrecep ya es fijo arriba */
  /*  Preparar salida: tamaño máximo posible considerando retardos */
  d = fmax((K / 2.0 + 3.0) * Lo + (Lo - 1.0), (K / 2.0 - 1.0) * Lo) +
      (Lo + 1151.0);
  i = Ry_gk->size[0] * Ry_gk->size[1];
  Ry_gk->size[0] = 1;
  b_loop_ub_tmp = (int)d;
  Ry_gk->size[1] = (int)d;
  emxEnsureCapacity_real_T(Ry_gk, i);
  setL1_data = Ry_gk->data;
  for (i = 0; i < b_loop_ub_tmp; i++) {
    setL1_data[i] = 0.0;
  }
  /*  BLOQUES 3 4: retardos, ponderaciones, combinación cruzada */
  i = (int)(K / 2.0);
  emxInit_real_T(&b_y0, 2);
  for (b_i = 0; b_i < i; b_i++) {
    double d1;
    double d2;
    int a_tmp;
    /*  Retardos  */
    d1 = (K / 2.0 + (4.0 - ((double)b_i + 1.0))) * Lo + (Lo - 1.0);
    d2 = (K / 2.0 - ((double)b_i + 1.0)) * Lo;
    /*  Inicializar  */
    i1 = x0->size[0] * x0->size[1];
    x0->size[0] = 1;
    x0->size[1] = (int)d;
    emxEnsureCapacity_real_T(x0, i1);
    x0_data = x0->data;
    i1 = b_y0->size[0] * b_y0->size[1];
    b_y0->size[0] = 1;
    b_y0->size[1] = (int)d;
    emxEnsureCapacity_real_T(b_y0, i1);
    setL0_data = b_y0->data;
    for (i1 = 0; i1 < b_loop_ub_tmp; i1++) {
      x0_data[i1] = 0.0;
      setL0_data[i1] = 0.0;
    }
    /*  Aplicar ponderaciones Hadamard */
    if (d1 + 1.0 > d1 + (Lo + 1151.0)) {
      i1 = 1;
    } else {
      i1 = (int)(d1 + 1.0);
    }
    a_tmp = H[((int)Code + (b_i << 2)) - 1];
    for (i2 = 0; i2 < loop_ub_tmp; i2++) {
      x0_data[(i1 + i2) - 1] = (double)a_tmp * C_data[(iv[b_i] + 2 * i2) - 1];
    }
    if (d2 + 1.0 > d2 + (Lo + 1151.0)) {
      i1 = 1;
    } else {
      i1 = (int)(d2 + 1.0);
    }
    for (i2 = 0; i2 < loop_ub_tmp; i2++) {
      setL0_data[(i1 + i2) - 1] =
          (double)a_tmp * S_data[(iv[b_i] + 2 * i2) - 1];
    }
    /*  Ventanas  */
    /*  Aplicación final */
    if (Ry_gk->size[1] == 1) {
      i1 = x0->size[1];
    } else {
      i1 = Ry_gk->size[1];
    }
    if ((Ry_gk->size[1] == x0->size[1]) && (i1 == b_y0->size[1])) {
      a_tmp = Ry_gk->size[1] - 1;
      i1 = Ry_gk->size[0] * Ry_gk->size[1];
      Ry_gk->size[0] = 1;
      emxEnsureCapacity_real_T(Ry_gk, i1);
      setL1_data = Ry_gk->data;
      for (i1 = 0; i1 <= a_tmp; i1++) {
        setL1_data[i1] = (setL1_data[i1] + x0_data[i1]) + setL0_data[i1];
      }
    } else {
      binary_expand_op(Ry_gk, x0, b_y0);
    }
  }
  emxFree_real_T(&b_y0);
  emxFree_real_T(&x0);
  emxFree_real_T(&S);
  emxFree_real_T(&C);
  /*  === Graficar señales internas === */
  /*  figure(2); plot(x0); grid on; title('Señal x0'); xlabel('Desplazamiento
   * (k)'); ylabel('Valor de Correlación'); */
  /*  figure(3); plot(y0); grid on; title('Señal y0'); xlabel('Desplazamiento
   * (k)'); ylabel('Valor de Correlación'); */
  /*  figure(4); plot(xT); grid on; title('Señal xT'); xlabel('Desplazamiento
   * (k)'); ylabel('Valor de Correlación'); */
  /*  figure(5); plot(yT); grid on; title('Señal yT'); xlabel('Desplazamiento
   * (k)'); ylabel('Valor de Correlación'); */
}

/* End of code generation (CorrLSmiracle_VELSC.c) */
