/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * correfMcss.c
 *
 * Code generation for function 'correfMcss'
 *
 */

/* Include files */
#include "correfMcss.h"
#include "CorrLSmiracle_VELSC_emxutil.h"
#include "CorrLSmiracle_VELSC_types.h"
#include "dec2bin.h"
#include "log2.h"
#include "mod.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void correfMcss(double M, const emxArray_real_T *x, double L, double W,
                emxArray_real_T *Rx)
{
  emxArray_char_T *bnr;
  emxArray_int8_T *coefn;
  emxArray_real_T *R;
  emxArray_real_T *Rdt;
  emxArray_real_T *b_W;
  emxArray_real_T *indn;
  const double *x_data;
  double *R_data;
  double *Rdt_data;
  double *Rx_data;
  double *W_data;
  double *indn_data;
  int b_i;
  int b_loop_ub_tmp;
  int i;
  int i1;
  int i2;
  int i3;
  int j;
  int k;
  int loop_ub;
  int loop_ub_tmp;
  char *bnr_data;
  signed char *coefn_data;
  x_data = x->data;
  /*  correfMcss realiza la correlacion eficiente con un conjunto de M
   * secuencias complementarias */
  /*  */
  /*  [Rx]=correfMcss(M,x,L,W); */
  /*  M:Numero de elementos del conjunto. */
  /*  x:Secuencia de entrada. */
  /*  L:longitud de las secuencias complementarias (M^N). */
  /*  W: Numero natural menor que L cuyo equivalente binario de n*N bits
   * representa  */
  /*  a las semillas que se utilizaron en la generacion de las secuencias
   * complementarias */
  /*  Rx: Salidas de la señal correlada con cada secuencia conjunto generada */
  /*  con la semilla W. */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  === CAMBIO MÍNIMO PARA MATLAB CODER === */
  /*  Asegura que Rx esté definido incluso si se entra en algún "return" de
   * validación */
  loop_ub_tmp = (int)M;
  i = Rx->size[0] * Rx->size[1];
  Rx->size[0] = (int)M;
  i1 = x->size[1];
  Rx->size[1] = x->size[1];
  emxEnsureCapacity_real_T(Rx, i);
  Rx_data = Rx->data;
  b_loop_ub_tmp = (int)M * x->size[1];
  for (i = 0; i < b_loop_ub_tmp; i++) {
    Rx_data[i] = 0.0;
  }
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  Chequeo de los argumentos de entrada % */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  emxInit_char_T(&bnr);
  emxInit_int8_T(&coefn);
  emxInit_real_T(&R, 2);
  emxInit_real_T(&Rdt, 2);
  emxInit_real_T(&indn, 2);
  emxInit_real_T(&b_W, 2);
  if (!(M < 2.0)) {
    double n;
    n = log(M);
    if (!(b_mod(n / 0.69314718055994529) != 0.0)) {
      double N;
      N = log(L) / n;
      if ((!(b_mod(N) != 0.0)) && (!(W >= L))) {
        double coef_tmp;
        int Rdt_tmp_tmp;
        int c_loop_ub_tmp;
        int d_loop_ub_tmp;
        signed char coef[4];
        /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         */
        /* Definicion de variables */
        /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         */
        n = b_log2(M);
        /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         */
        /*  Conversión de W a vector binario ANTES de indexarlo */
        /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         */
        dec2bin(W, n * N, bnr);
        bnr_data = bnr->data;
        i = b_W->size[0] * b_W->size[1];
        b_W->size[0] = 1;
        loop_ub = bnr->size[1];
        b_W->size[1] = bnr->size[1];
        emxEnsureCapacity_real_T(b_W, i);
        W_data = b_W->data;
        for (i = 0; i < loop_ub; i++) {
          char c;
          c = bnr_data[i];
          W_data[i] = (c == '1') - (c != '1');
        }
        /*  transforma los ceros en -1s */
        /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         */
        /* Generacion de coeficientes para  */
        /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         */
        /* Generacion de coeficientes para cada etapa */
        /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         */
        i = coefn->size[0] * coefn->size[1] * coefn->size[2];
        coefn->size[0] = (int)M;
        coefn->size[1] = (int)M;
        c_loop_ub_tmp = (int)N;
        coefn->size[2] = (int)N;
        emxEnsureCapacity_int8_T(coefn, i);
        coefn_data = coefn->data;
        d_loop_ub_tmp = (int)M * (int)M * (int)N;
        for (i = 0; i < d_loop_ub_tmp; i++) {
          coefn_data[i] = 0;
        }
        coef[0] = 1;
        coef[2] = (signed char)W_data[0];
        coef[1] = 1;
        coef[3] = (signed char)-(int)W_data[0];
        for (i = 0; i < loop_ub_tmp; i++) {
          for (i2 = 0; i2 < loop_ub_tmp; i2++) {
            coefn_data[i2 + coefn->size[0] * i] = coef[i2 + (int)M * i];
          }
        }
        i = (int)(N - 1.0);
        if ((int)(N - 1.0) - 1 >= 0) {
          coef[0] = 1;
          coef[1] = 1;
        }
        for (j = 0; j < i; j++) {
          coef_tmp = W_data[(int)((n * ((double)j + 2.0) - n) + 1.0) - 1];
          coef[2] = (signed char)coef_tmp;
          coef[3] = (signed char)-(int)coef_tmp;
          for (i2 = 0; i2 < loop_ub_tmp; i2++) {
            for (i3 = 0; i3 < loop_ub_tmp; i3++) {
              coefn_data[(i3 + coefn->size[0] * i2) +
                         coefn->size[0] * coefn->size[1] * (j + 1)] =
                  coef[i3 + (int)M * i2];
            }
          }
        }
        /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         */
        /* Inicializacion de variables */
        /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         */
        /*  Valores acumulados (Buffers de memoria) */
        d_loop_ub_tmp = (int)(M - 1.0);
        i2 = R->size[0] * R->size[1];
        R->size[0] = (int)(M - 1.0);
        R->size[1] = (int)L;
        emxEnsureCapacity_real_T(R, i2);
        R_data = R->data;
        loop_ub = (int)(M - 1.0) * (int)L;
        for (i2 = 0; i2 < loop_ub; i2++) {
          R_data[i2] = 0.0;
        }
        /*  buffer fijo */
        /*  Inicializacion de las salidas */
        i2 = Rx->size[0] * Rx->size[1];
        Rx->size[0] = (int)M;
        Rx->size[1] = x->size[1];
        emxEnsureCapacity_real_T(Rx, i2);
        Rx_data = Rx->data;
        for (i2 = 0; i2 < b_loop_ub_tmp; i2++) {
          Rx_data[i2] = 0.0;
        }
        /*  salida fija: M x length(x) */
        /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         */
        /* Algoritmo principal */
        /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         */
        /* Itera la señal de entrada  */
        b_loop_ub_tmp = (int)M * (int)N;
        Rdt_tmp_tmp = (int)M - 1;
        for (b_i = 0; b_i < i1; b_i++) {
          double retardo;
          i2 = b_W->size[0] * b_W->size[1];
          b_W->size[0] = 1;
          b_W->size[1] = (int)(M - 1.0);
          emxEnsureCapacity_real_T(b_W, i2);
          W_data = b_W->data;
          for (i2 = 0; i2 < d_loop_ub_tmp; i2++) {
            W_data[i2] = 0.0;
          }
          i2 = Rdt->size[0] * Rdt->size[1];
          Rdt->size[0] = (int)M;
          Rdt->size[1] = (int)N;
          emxEnsureCapacity_real_T(Rdt, i2);
          Rdt_data = Rdt->data;
          for (i2 = 0; i2 < b_loop_ub_tmp; i2++) {
            Rdt_data[i2] = 0.0;
          }
          /* Etapa 1 */
          retardo = L / M;
          /*  primera etapa de retardos */
          for (j = 0; j < d_loop_ub_tmp; j++) {
            W_data[j] = retardo * (M - ((double)j + 1.0)) + 1.0;
          }
          /*  Actualiza retardos (shift register) */
          for (k = 0; k < d_loop_ub_tmp; k++) {
            if (R->size[1] < 2) {
              i2 = 0;
              i3 = 0;
            } else {
              i2 = 1;
              i3 = (int)L;
            }
            loop_ub = i3 - i2;
            i3 = indn->size[0] * indn->size[1];
            indn->size[0] = 1;
            indn->size[1] = loop_ub;
            emxEnsureCapacity_real_T(indn, i3);
            indn_data = indn->data;
            for (i3 = 0; i3 < loop_ub; i3++) {
              indn_data[i3] = R_data[k + R->size[0] * i3];
            }
            for (i3 = 0; i3 < loop_ub; i3++) {
              R_data[k + R->size[0] * (i2 + i3)] = indn_data[i3];
            }
            R_data[k] = x_data[b_i];
          }
          /*  Correlación jerárquica */
          n = R_data[(int)M - 2];
          Rdt_data[(int)M - 1] =
              (double)coefn_data[((int)M + coefn->size[0] * ((int)M - 1)) - 1] *
              n;
          for (k = 0; k <= loop_ub_tmp - 2; k++) {
            Rdt_data[(int)M - 1] +=
                (double)coefn_data[((int)M + coefn->size[0] * k) - 1] *
                R_data[k + R->size[0] * ((int)W_data[k] - 1)];
          }
          for (k = 0; k <= loop_ub_tmp - 2; k++) {
            for (i2 = 0; i2 < c_loop_ub_tmp; i2++) {
              Rdt_data[k + Rdt->size[0] * i2] =
                  (double)coefn_data[k + coefn->size[0] * ((int)M - 1)] * n;
            }
            for (loop_ub = 0; loop_ub <= loop_ub_tmp - 2; loop_ub++) {
              coef_tmp =
                  (double)coefn_data[k + coefn->size[0] * loop_ub] *
                  R_data[loop_ub + R->size[0] * ((int)W_data[loop_ub] - 1)];
              for (i2 = 0; i2 < c_loop_ub_tmp; i2++) {
                Rdt_data[k + Rdt->size[0] * i2] += coef_tmp;
              }
            }
          }
          /*  ACTUALIZACIÓN DE RETARDOS (ESCALAR) */
          for (k = 0; k <= loop_ub_tmp - 2; k++) {
            R_data[k + R->size[0] * ((int)W_data[k] - 1)] = Rdt_data[k];
          }
          /* Etapas 2 hasta N */
          for (j = 0; j < i; j++) {
            retardo /= M;
            i2 = indn->size[0] * indn->size[1];
            indn->size[0] = 1;
            indn->size[1] = (int)M - 1;
            emxEnsureCapacity_real_T(indn, i2);
            indn_data = indn->data;
            for (k = 0; k < Rdt_tmp_tmp; k++) {
              indn_data[k] = W_data[k] + retardo * (M - ((double)k + 1.0));
            }
            i2 = b_W->size[0] * b_W->size[1];
            b_W->size[0] = 1;
            b_W->size[1] = (int)M - 1;
            emxEnsureCapacity_real_T(b_W, i2);
            W_data = b_W->data;
            for (i2 = 0; i2 < Rdt_tmp_tmp; i2++) {
              W_data[i2] = indn_data[i2];
            }
            Rdt_data[((int)M + Rdt->size[0] * (j + 1)) - 1] =
                (double)coefn_data[(((int)M + coefn->size[0] * ((int)M - 1)) +
                                    coefn->size[0] * coefn->size[1] * (j + 1)) -
                                   1] *
                Rdt_data[((int)M + Rdt->size[0] * j) - 1];
            for (k = 0; k <= loop_ub_tmp - 2; k++) {
              Rdt_data[((int)M + Rdt->size[0] * (j + 1)) - 1] +=
                  (double)
                      coefn_data[(((int)M + coefn->size[0] * k) +
                                  coefn->size[0] * coefn->size[1] * (j + 1)) -
                                 1] *
                  R_data[k + R->size[0] * ((int)indn_data[k] - 1)];
            }
            for (k = 0; k <= loop_ub_tmp - 2; k++) {
              n = Rdt_data[((int)M + Rdt->size[0] * j) - 1];
              for (i2 = 0; i2 < c_loop_ub_tmp; i2++) {
                Rdt_data[k + Rdt->size[0] * i2] =
                    (double)
                        coefn_data[(k + coefn->size[0] * ((int)M - 1)) +
                                   coefn->size[0] * coefn->size[1] * (j + 1)] *
                    n;
              }
              for (loop_ub = 0; loop_ub <= loop_ub_tmp - 2; loop_ub++) {
                coef_tmp =
                    (double)
                        coefn_data[(k + coefn->size[0] * loop_ub) +
                                   coefn->size[0] * coefn->size[1] * (j + 1)] *
                    R_data[loop_ub +
                           R->size[0] * ((int)indn_data[loop_ub] - 1)];
                for (i2 = 0; i2 < c_loop_ub_tmp; i2++) {
                  Rdt_data[k + Rdt->size[0] * i2] += coef_tmp;
                }
              }
            }
            for (k = 0; k <= loop_ub_tmp - 2; k++) {
              R_data[k + R->size[0] * ((int)indn_data[k] - 1)] =
                  Rdt_data[k + Rdt->size[0] * (j + 1)];
            }
          }
          /*  Resultado final */
          for (j = 0; j < loop_ub_tmp; j++) {
            if (j + 1 == (int)M) {
              Rx_data[j + Rx->size[0] * b_i] =
                  Rdt_data[((int)M + Rdt->size[0] * ((int)N - 1)) - 1];
            } else {
              Rx_data[j + Rx->size[0] * b_i] =
                  R_data[j + R->size[0] * ((int)W_data[j] - 1)];
            }
          }
        }
      }
    }
  }
  emxFree_real_T(&b_W);
  emxFree_real_T(&indn);
  emxFree_real_T(&Rdt);
  emxFree_real_T(&R);
  emxFree_int8_T(&coefn);
  emxFree_char_T(&bnr);
}

/* End of code generation (correfMcss.c) */
