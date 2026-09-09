/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * dec2bin.c
 *
 * Code generation for function 'dec2bin'
 *
 */

/* Include files */
#include "dec2bin.h"
#include "CorrLSmiracle_VELSC_emxutil.h"
#include "CorrLSmiracle_VELSC_types.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
void dec2bin(double d, double numBits, emxArray_char_T *s)
{
  double y;
  int idx;
  int nbits;
  int p;
  int pmax;
  int pmin;
  int pow2p;
  char sfull[64];
  char *s_data;
  y = numBits;
  if (fabs(numBits) < 4.503599627370496E+15) {
    if (numBits >= 0.5) {
      y = floor(numBits + 0.5);
    } else if (numBits > -0.5) {
      y = numBits * 0.0;
    } else {
      y = ceil(numBits - 0.5);
    }
  }
  if ((int)y <= 1) {
    nbits = 1;
  } else {
    nbits = (int)y;
  }
  if (d < 0.0) {
    char b[64];
    boolean_T carry;
    boolean_T exitg1;
    y = -d;
    for (p = 0; p < 64; p++) {
      b[p] = '0';
    }
    idx = 64;
    exitg1 = false;
    while ((!exitg1) && (idx > 0)) {
      double olddi;
      olddi = y;
      y /= 2.0;
      y = floor(y);
      if (2.0 * y < olddi) {
        b[idx - 1] = '1';
      }
      if (!(y > 0.0)) {
        exitg1 = true;
      } else {
        idx--;
      }
    }
    for (pow2p = 0; pow2p < 64; pow2p++) {
      sfull[pow2p] = '1';
      if (b[pow2p] == '1') {
        sfull[pow2p] = '0';
      }
    }
    carry = true;
    for (pow2p = 63; pow2p >= 0; pow2p--) {
      if (carry) {
        if (sfull[pow2p] == '1') {
          sfull[pow2p] = '0';
        } else {
          sfull[pow2p] = '1';
          carry = false;
        }
      }
    }
    idx = 0;
    pow2p = 0;
    exitg1 = false;
    while ((!exitg1) && (pow2p < 64)) {
      if (sfull[pow2p] == '0') {
        idx = pow2p + 1;
        exitg1 = true;
      } else {
        pow2p++;
      }
    }
    if (idx == 0) {
      pmin = 55;
      pmax = 8;
    } else {
      idx--;
      if (65 - idx <= 4) {
        idx = 4;
      } else {
        idx = 65 - idx;
      }
      pmax = 31;
      pmin = 0;
      exitg1 = false;
      while ((!exitg1) && (pmax - pmin > 1)) {
        p = (pmin + pmax) >> 1;
        pow2p = 1 << p;
        if (pow2p == idx) {
          pmax = p;
          exitg1 = true;
        } else if (pow2p > idx) {
          pmax = p;
        } else {
          pmin = p;
        }
      }
      idx = 1 << pmax;
      if (idx >= 64) {
        pmax = 64;
      } else {
        pmax = idx;
      }
      pmin = 63 - pmax;
    }
  } else {
    boolean_T exitg1;
    y = d;
    for (p = 0; p < 64; p++) {
      sfull[p] = '0';
    }
    idx = 64;
    exitg1 = false;
    while ((!exitg1) && (idx > 0)) {
      double olddi;
      olddi = y;
      y /= 2.0;
      y = floor(y);
      if (2.0 * y < olddi) {
        sfull[idx - 1] = '1';
      }
      if (!(y > 0.0)) {
        exitg1 = true;
      } else {
        idx--;
      }
    }
    idx = 0;
    pow2p = 0;
    exitg1 = false;
    while ((!exitg1) && (pow2p < 64)) {
      if (sfull[pow2p] == '1') {
        idx = pow2p + 1;
        exitg1 = true;
      } else {
        pow2p++;
      }
    }
    if (idx == 0) {
      idx = 64;
    }
    pmin = idx - 2;
    pmax = 65 - idx;
  }
  if (nbits >= pmax) {
    p = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = nbits;
    emxEnsureCapacity_char_T(s, p);
    s_data = s->data;
  } else {
    p = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = pmax;
    emxEnsureCapacity_char_T(s, p);
    s_data = s->data;
  }
  if (pmax < nbits) {
    char padval;
    idx = (nbits - pmax) - 1;
    if (d < 0.0) {
      padval = '1';
    } else {
      padval = '0';
    }
    for (pow2p = 0; pow2p <= idx; pow2p++) {
      s_data[pow2p] = padval;
    }
  } else {
    idx = -1;
  }
  p = (unsigned char)pmax;
  for (pow2p = 0; pow2p < p; pow2p++) {
    s_data[(idx + pow2p) + 1] = sfull[(pmin + pow2p) + 1];
  }
}

/* End of code generation (dec2bin.c) */
