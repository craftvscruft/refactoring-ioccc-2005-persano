#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const int H = 20;

typedef double double_4[4];

double o[3], d[3], v[3];
double (*q)[20][6];
double C[80][20][6];
double sr = 15;
double br = 25;
double A;
double s;
double t;
const int e = 256;

int D[64];
char B[256];

void T(int *p, int k, int n) {
    int g;

    *p = n;

    for (g = 0; k && g < 4; g++) {
        T(p + (9 & (g * 4 | g)) * (1 << k - 1), k - 1,
          n + (1 << 6 - 2 * k) * (4 - g & 3));
    }
}

void O(int c) {
    static int S;
    B[++S] = c,

            (S > e - 2 || c == 129) && (*B = S, S = fwrite(B, S + 1, 1, stdout) - 1);
}

void N(double *v) {
    double w = sqrt(*v * *v + 1[v] * v[1] + v[2] * 2[v]);
    for (int k = 0; k < 3; k++) {
        v[k] = v[k] / w;
    }
}

void K(int j) {
    double n[3], p[3];
    double *b = (*q)[j];
    double w = cos(A);
    double x = sin(A);
    int k;
    for (k = 0; k < 3; k++) {
        p[k] = br * o[k] +
               sr * (n[k] = cos(t) * k[v] + sin(t) * (v[(k + 1) % 3] * d[(k + 2) % 3] -
                                                  v[(k + 2) % 3] * d[(k + 1) % 3]));
    }
    *b++ = p[1] + e / 2, *b++ = -*p * x + p[2] * w - e * e,
    *b++ = (*p * w + p[2] * x) + e / 2;
    for (k = 0; k < 3; k++) {
        b[k] = fabs(*n * !!k + n[1] * (k < 2) + n[2]) / sqrt(2 + !(k - 1));
    }
}
double getD(int s) {
    return pow(25 / log(36), 4) / s / 377;
}



void m_loop_K(double *a, int s) {
    double d = getD(s);
    for (int i = *a = 0; i < s; i++, *a += d) {
        K(i);
    }
}

void G(int i, int P, int Q) {
    double w = -sin(P * s);
    double x = cos(P * s);
    double y = sin(Q * s);
    double z = cos(Q * s);
    int k;
    for (k = 0; k < 3; k++) {
        v[k] = k[o] = (k < 2) * (x + 3) * (k ? y : z) - !(k - 2) * w;
    }
    N(v);
    for (k = 0; k < 3; k++) {
        d[k] =
                P * (k ? k - 1 ? x : w * y : w * z) + Q * (k ? k - 1 ? 0 : *o : -o[1]);
    }
    N(d), w = *d * *v + d[1] * v[1] + d[2] * v[2];
    for (k = 0; k < 3; k++) {
        v[k] = k[v] - w * d[k];
    }
    N(v), q = &C[i], m_loop_K(&t, H);
}

void W(char *s) { *s && (W(s + 1), putchar(*s - 98 * (*s > 97))); }

void m_loop_G(double *a, int s, int P, int Q) {
    double d = getD(s);
    for (int i = *a = 0; i < s; i++, *a += d) {
        G(i, P, Q);
    }
}


void E(int z, int P, int Q) {
    double *r;
    double *a;
    double *b;
    double_4 *q, l, x, d, I[256];
    int i;
    int j;
    int m;
    int c;
    int y;
    int w;
    int h;
    int Y = sizeof l;
    fputc(46, stderr), m_loop_G(&s, 4 * H, P, Q), W("ibcbcbbbbb,");
    int g = 130;
    O(e / 2);
    for (y = 0; y < e; y++) {
        for (i = e; i;) {
            *(I[--i]) = 0;
        }
        int k;
        for (; i < 4 * H * H; i++) {
            for (*l = k = 0; k - 5; k++, a = b) {
                if (b = C[(i / H + ((k ^ k / 2) & 1)) % (4 * H)][(i + (k / 2 & 1)) % H],
                        k && y < *a ^ y < *b) {
                    for (h = 0; h < 4; h++) {
                        h[x] = a[h + 2] + (b[h + 2] - a[h + 2]) * (y - *a) / (*b - *a),
                                h ? *l && (d[h] = (l[h] - h[x]) / w)
                                  : (w = 1 + fabs(*l - (*x = (int) *x)));
                    }
                    for (q = I + (int) *x; *l && w--; q += 2 * (*x < *l) - 1) {
                        for (**q > C[i / H][i % H][1] &&
                             (memcpy(*q, x, Y), **q = C[i / H][i % H][1]),
                                     h = 1;
                             h < 4; h++) {
                            x[h] += h[d];
                        }
                    }
                    memcpy(l, x, Y);
                }
            }
        }
        for (i = 0; i < e;
             i++, O(*r ? c : 127), g = g < e - 1 ? g + 1 : (O(e / 2), 130)) {
            for (r = &I[i][3], c = 0, j = 2; j + 1; j--, r--) {
                k = 3 | !(j - 1) * 4,
                m = *r * k +
                    (63 * (*r * k - (int) (*r * k)) > D[(y * 8 & 56) + (i & 7)]),
                c <<= k / 4 + 2, c |= (m | -(m > k)) & k;
            }
        }
    }
    O(129);
    putchar(0);
}

void m_loopE(double *a, int s, int P, int Q) {
    double d = getD(s);
    for (int i = *a = 0; i < s; i++, *a += d) {
        E(i, P, Q);
    }
}





int main(int c, char **v) {
    if (c < 3) {
        return fprintf(stderr, "Usage: %s P Q [F]\n", *v), 1;
    }

    int P = atoi(v[1]);
    int Q = atoi(v[2]);
    T(D, 3, 0);

    W("bb\346cbcba98FIG");

    for (int k=0; k < 384; k++) {
        putchar(k / 3 << ("gdb"[k % 3] - 97) | 31 | (k % 3 != 1) << 5);
    }

    W("bbbce0.2EPACSTEN\x0b\xff!"),

            m_loopE(&A, c > 3 ? atoi(v[3]) : 40, P, Q),

            putchar(59);

    return 0;
}
