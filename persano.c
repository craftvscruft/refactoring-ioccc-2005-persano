#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const int H = 20;

typedef double R, V[3], U[4], (*J)(R);

V o, d, v;
R (*q)[H][6], C[4 * H][H][6], sr = 15, br = 25, A, s, t;
int P;
int Q;
int g;
int S;
int k;
int e = 256;

int (*u)(int) = putchar;

int D[64];
J f = cos, Y = sin;
char B[256];

void T(int *p, int k, int n) {
  int g;

  *p = n;

  for (g = 0; k && g < 4; g++)
    T(p + (9 & (g * 4 | g)) * (1 << k - 1), k - 1,
      n + (1 << 6 - 2 * k) * (4 - g & 3));
}

void O(int c) {
  B[++S] = c,

  (S > e - 2 || c == 129) && (*B = S, S = fwrite(B, S + 1, 1, stdout) - 1);
}

void N(R *v) {
  R w = sqrt(*v * *v + 1 [v] * v[1] + v[2] * 2 [v]);
  for (k = 0; k < 3; k++)
    v[k] = v[k] / w;
}
void K(int j) {
  V n, p;
  R *b = (*q)[j], w = f(A), x = Y(A);
  for (k = 0; k < 3; k++)
    p[k] = br * o[k] +
           sr * (n[k] = f(t) * k[v] + Y(t) * (v[(k + 1) % 3] * d[(k + 2) % 3] -
                                              v[(k + 2) % 3] * d[(k + 1) % 3]));
  *b++ = p[1] + e / 2, *b++ = -*p * x + p[2] * w - e * e,
  *b++ = (*p * w + p[2] * x) + e / 2;
  for (k = 0; k < 3; k++)
    b[k] = fabs(*n * !!k + n[1] * (k < 2) + n[2]) / sqrt(2 + !(k - 1));
}
void M(R *a, int s, void (*m)(int)) {
  int i;
  R d = pow(25 / log(36), 4) / s / 377;
  for (*a = i = 0; i < s; i++, *a += d)
    m(i);
}
void G(int i) {
  R w = -Y(P * s), x = f(P * s), y = Y(Q * s), z = f(Q * s);
  for (k = 0; k < 3; k++)
    v[k] = k[o] = (k < 2) * (x + 3) * (k ? y : z) - !(k - 2) * w;
  N(v);
  for (k = 0; k < 3; k++)
    d[k] =
        P * (k ? k - 1 ? x : w * y : w * z) + Q * (k ? k - 1 ? 0 : *o : -o[1]);
  N(d), w = *d * *v + d[1] * v[1] + d[2] * v[2];
  for (k = 0; k < 3; k++)
    v[k] = k[v] - w * d[k];
  N(v), q = &C[i], M(&t, H, K);
}
void W(char *s) { *s && (W(s + 1), u(*s - 98 * (*s > 97))); }
void E(int z) {
  R *r, *a, *b;
  U *q, l, x, d, I[256];
  int i;
  int j;
  int m;
  int c;
  int y;
  int w;
  int h;
  int Y = sizeof l;
  fputc(46, stderr), M(&s, 4 * H, G), W("ibcbcbbbbb,"), g = 130, O(e / 2);
  for (y = 0; y < e; y++) {
    for (i = e; i;)
      *(I[--i]) = 0;
    for (; i < 4 * H * H; i++)
      for (*l = k = 0; k - 5; k++, a = b)
        if (b = C[(i / H + ((k ^ k / 2) & 1)) % (4 * H)][(i + (k / 2 & 1)) % H],
            k && y < *a ^ y < *b) {
          for (h = 0; h < 4; h++)
            h[x] = a[h + 2] + (b[h + 2] - a[h + 2]) * (y - *a) / (*b - *a),
            h ? *l && (d[h] = (l[h] - h[x]) / w)
              : (w = 1 + fabs(*l - (*x = (int)*x)));
          for (q = I + (int)*x; *l && w--; q += 2 * (*x < *l) - 1)
            for (**q > C[i / H][i % H][1] &&
                     (memcpy(*q, x, Y), **q = C[i / H][i % H][1]),
                 h = 1;
                 h < 4; h++)
              x[h] += h[d];
          memcpy(l, x, Y);
        }
    for (i = 0; i < e;
         i++, O(*r ? c : 127), g = g < e - 1 ? g + 1 : (O(e / 2), 130))
      for (r = &I[i][3], c = 0, j = 2; j + 1; j--, r--)
        k = 3 | !(j - 1) * 4,
        m = *r * k +
            (63 * (*r * k - (int)(*r * k)) > D[(y * 8 & 56) + (i & 7)]),
        c <<= k / 4 + 2, c |= (m | -(m > k)) & k;
  }
  O(129);
  u(0);
}
int
    
    
	
main(int c, char **v)
{
  if (c < 3)
    return fprintf(stderr, "Usage: %s P Q [F]\n", *v), 1;

  P = atoi(v[1]), Q = atoi(v[2]),

  T(D, 3, 0);

  W("bb\346cbcba98FIG");

  for (; k < 384; k++)
    u(k / 3 << ("gdb"[k % 3] - 97) | 31 | (k % 3 != 1) << 5);

  W("bbbce0.2EPACSTEN\x0b\xff!"),

      M(&A, c > 3 ? atoi(v[3]) : 40, E),

      u(59);

  return 0;
}

/* vi:set ts=8: */
