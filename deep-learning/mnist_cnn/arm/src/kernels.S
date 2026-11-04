// AArch64 NEON kernels. Each processes 4 floats per loop, then a scalar tail.
// Only v0-v7 and v16-v31 are used, because v8-v15 are callee-saved.
// Note: on Darwin arm64 a semicolon starts a comment, so one instruction per line.

.text

.macro FUNC name
    .globl _\name
    .p2align 2
_\name:
.endm

// void k_saxpy(long n, float a, const float *x, float *y)   y += a * x
FUNC k_saxpy
    dup     v0.4s, v0.s[0]
1:  cmp     x0, #4
    b.lt    2f
    ldr     q1, [x1], #16
    ldr     q2, [x2]
    fmla    v2.4s, v1.4s, v0.4s
    str     q2, [x2], #16
    sub     x0, x0, #4
    b       1b
2:  cbz     x0, 3f
    ldr     s1, [x1], #4
    ldr     s2, [x2]
    fmadd   s2, s1, s0, s2
    str     s2, [x2], #4
    sub     x0, x0, #1
    b       2b
3:  ret

// float k_sdot(long n, const float *x, const float *y)
FUNC k_sdot
    movi    v0.4s, #0
1:  cmp     x0, #4
    b.lt    2f
    ldr     q1, [x1], #16
    ldr     q2, [x2], #16
    fmla    v0.4s, v1.4s, v2.4s
    sub     x0, x0, #4
    b       1b
2:  faddp   v0.4s, v0.4s, v0.4s
    faddp   s0, v0.2s
3:  cbz     x0, 4f
    ldr     s1, [x1], #4
    ldr     s2, [x2], #4
    fmadd   s0, s1, s2, s0
    sub     x0, x0, #1
    b       3b
4:  ret

// float k_ssum(long n, const float *x)
FUNC k_ssum
    movi    v0.4s, #0
1:  cmp     x0, #4
    b.lt    2f
    ldr     q1, [x1], #16
    fadd    v0.4s, v0.4s, v1.4s
    sub     x0, x0, #4
    b       1b
2:  faddp   v0.4s, v0.4s, v0.4s
    faddp   s0, v0.2s
3:  cbz     x0, 4f
    ldr     s1, [x1], #4
    fadd    s0, s0, s1
    sub     x0, x0, #1
    b       3b
4:  ret

// void k_scale_shift(long n, float a, float b, const float *src, float *dst)   dst = a * src + b
FUNC k_scale_shift
    dup     v0.4s, v0.s[0]
    dup     v1.4s, v1.s[0]
1:  cmp     x0, #4
    b.lt    2f
    ldr     q2, [x1], #16
    mov     v3.16b, v1.16b
    fmla    v3.4s, v2.4s, v0.4s
    str     q3, [x2], #16
    sub     x0, x0, #4
    b       1b
2:  cbz     x0, 3f
    ldr     s2, [x1], #4
    fmadd   s3, s2, s0, s1
    str     s3, [x2], #4
    sub     x0, x0, #1
    b       2b
3:  ret

// void k_relu(long n, const float *src, float *dst)
FUNC k_relu
    movi    v0.4s, #0
1:  cmp     x0, #4
    b.lt    2f
    ldr     q1, [x1], #16
    fmax    v1.4s, v1.4s, v0.4s
    str     q1, [x2], #16
    sub     x0, x0, #4
    b       1b
2:  cbz     x0, 3f
    ldr     s1, [x1], #4
    fmax    s1, s1, s0
    str     s1, [x2], #4
    sub     x0, x0, #1
    b       2b
3:  ret

// void k_relu_back(long n, const float *gy, const float *x, float *gx)   gx = gy where x > 0, else 0
FUNC k_relu_back
1:  cmp     x0, #4
    b.lt    2f
    ldr     q1, [x1], #16
    ldr     q2, [x2], #16
    fcmgt   v2.4s, v2.4s, #0.0
    and     v1.16b, v1.16b, v2.16b
    str     q1, [x3], #16
    sub     x0, x0, #4
    b       1b
2:  cbz     x0, 3f
    ldr     s1, [x1], #4
    ldr     s2, [x2], #4
    fcmgt   s2, s2, #0.0
    and     v1.8b, v1.8b, v2.8b
    str     s1, [x3], #4
    sub     x0, x0, #1
    b       2b
3:  ret

// void k_conv3_row(long n, float *dst, const float *r0, const float *r1, const float *r2, const float *w)
// dst[i] += sum over k,l of w[3k+l] * rk[i-1+l]. Rows must be readable one float either side.
.macro TAPS4 row, wa, wb, wc
    ldur    q1, [\row, #-4]
    fmla    v0.4s, v1.4s, \wa\().4s
    ldr     q1, [\row]
    fmla    v0.4s, v1.4s, \wb\().4s
    ldur    q1, [\row, #4]
    fmla    v0.4s, v1.4s, \wc\().4s
    add     \row, \row, #16
.endm
.macro TAPS1 row, wa, wb, wc
    ldur    s1, [\row, #-4]
    fmadd   s0, s1, \wa, s0
    ldr     s1, [\row]
    fmadd   s0, s1, \wb, s0
    ldur    s1, [\row, #4]
    fmadd   s0, s1, \wc, s0
    add     \row, \row, #4
.endm
FUNC k_conv3_row
    ld1r    {v16.4s}, [x5], #4
    ld1r    {v17.4s}, [x5], #4
    ld1r    {v18.4s}, [x5], #4
    ld1r    {v19.4s}, [x5], #4
    ld1r    {v20.4s}, [x5], #4
    ld1r    {v21.4s}, [x5], #4
    ld1r    {v22.4s}, [x5], #4
    ld1r    {v23.4s}, [x5], #4
    ld1r    {v24.4s}, [x5], #4
1:  cmp     x0, #4
    b.lt    2f
    ldr     q0, [x1]
    TAPS4   x2, v16, v17, v18
    TAPS4   x3, v19, v20, v21
    TAPS4   x4, v22, v23, v24
    str     q0, [x1], #16
    sub     x0, x0, #4
    b       1b
2:  cbz     x0, 3f
    ldr     s0, [x1]
    TAPS1   x2, s16, s17, s18
    TAPS1   x3, s19, s20, s21
    TAPS1   x4, s22, s23, s24
    str     s0, [x1], #4
    sub     x0, x0, #1
    b       2b
3:  ret

// void k_conv3_grad_row(long n, float *wg, const float *g, const float *r0, const float *r1, const float *r2)
// wg[3k+l] += sum over i of g[i] * rk[i-1+l]
.macro GRAD4 row, a, b, c
    ldur    q1, [\row, #-4]
    fmla    \a\().4s, v1.4s, v0.4s
    ldr     q1, [\row]
    fmla    \b\().4s, v1.4s, v0.4s
    ldur    q1, [\row, #4]
    fmla    \c\().4s, v1.4s, v0.4s
    add     \row, \row, #16
.endm
.macro GRAD1 row, a, b, c
    ldur    s1, [\row, #-4]
    fmadd   \a, s1, s0, \a
    ldr     s1, [\row]
    fmadd   \b, s1, s0, \b
    ldur    s1, [\row, #4]
    fmadd   \c, s1, s0, \c
    add     \row, \row, #4
.endm
.macro HSUM v, s
    faddp   \v\().4s, \v\().4s, \v\().4s
    faddp   \s, \v\().2s
.endm
.macro ADDTO s, off
    ldr     s2, [x1, #\off]
    fadd    s2, s2, \s
    str     s2, [x1, #\off]
.endm
FUNC k_conv3_grad_row
    movi    v16.4s, #0
    movi    v17.4s, #0
    movi    v18.4s, #0
    movi    v19.4s, #0
    movi    v20.4s, #0
    movi    v21.4s, #0
    movi    v22.4s, #0
    movi    v23.4s, #0
    movi    v24.4s, #0
1:  cmp     x0, #4
    b.lt    2f
    ldr     q0, [x2], #16
    GRAD4   x3, v16, v17, v18
    GRAD4   x4, v19, v20, v21
    GRAD4   x5, v22, v23, v24
    sub     x0, x0, #4
    b       1b
2:  HSUM    v16, s16
    HSUM    v17, s17
    HSUM    v18, s18
    HSUM    v19, s19
    HSUM    v20, s20
    HSUM    v21, s21
    HSUM    v22, s22
    HSUM    v23, s23
    HSUM    v24, s24
3:  cbz     x0, 4f
    ldr     s0, [x2], #4
    GRAD1   x3, s16, s17, s18
    GRAD1   x4, s19, s20, s21
    GRAD1   x5, s22, s23, s24
    sub     x0, x0, #1
    b       3b
4:  ADDTO   s16, 0
    ADDTO   s17, 4
    ADDTO   s18, 8
    ADDTO   s19, 12
    ADDTO   s20, 16
    ADDTO   s21, 20
    ADDTO   s22, 24
    ADDTO   s23, 28
    ADDTO   s24, 32
    ret

// void k_maxpool2_row(long n_out, const float *r0, const float *r1, float *dst, int32_t *idx)
// idx is 0..3 for r0[2i], r0[2i+1], r1[2i], r1[2i+1]. The first maximum wins, as in PyTorch.
FUNC k_maxpool2_row
    movi    v16.4s, #1
    movi    v17.4s, #2
    movi    v18.4s, #3
1:  cmp     x0, #4
    b.lt    2f
    ld2     {v0.4s, v1.4s}, [x1], #32
    ld2     {v2.4s, v3.4s}, [x2], #32
    fcmgt   v4.4s, v1.4s, v0.4s
    and     v6.16b, v4.16b, v16.16b
    bsl     v4.16b, v1.16b, v0.16b
    fcmgt   v5.4s, v2.4s, v4.4s
    bit     v4.16b, v2.16b, v5.16b
    bit     v6.16b, v17.16b, v5.16b
    fcmgt   v5.4s, v3.4s, v4.4s
    bit     v4.16b, v3.16b, v5.16b
    bit     v6.16b, v18.16b, v5.16b
    str     q4, [x3], #16
    str     q6, [x4], #16
    sub     x0, x0, #4
    b       1b
2:  cbz     x0, 3f
    ldp     s0, s1, [x1], #8
    ldp     s2, s3, [x2], #8
    mov     w5, #0
    mov     w6, #1
    fcmp    s1, s0
    fcsel   s0, s1, s0, gt
    csel    w5, w6, w5, gt
    mov     w6, #2
    fcmp    s2, s0
    fcsel   s0, s2, s0, gt
    csel    w5, w6, w5, gt
    mov     w6, #3
    fcmp    s3, s0
    fcsel   s0, s3, s0, gt
    csel    w5, w6, w5, gt
    str     s0, [x3], #4
    str     w5, [x4], #4
    sub     x0, x0, #1
    b       2b
3:  ret

// void k_adamw(long n, float *w, const float *g, float *m, float *v, float b1, float b2,
//              float eps, float decay, float step, float root_inv)
// w *= decay, m = b1 m + (1-b1) g, v = b2 v + (1-b2) g^2, w -= step * m / (sqrt(v) * root_inv + eps)
FUNC k_adamw
    fmov    s6, #1.0
    fsub    s7, s6, s0
    fsub    s6, s6, s1
    dup     v16.4s, v0.s[0]
    dup     v17.4s, v1.s[0]
    dup     v18.4s, v2.s[0]
    dup     v19.4s, v3.s[0]
    dup     v20.4s, v4.s[0]
    dup     v21.4s, v5.s[0]
    dup     v22.4s, v7.s[0]
    dup     v23.4s, v6.s[0]
1:  cmp     x0, #4
    b.lt    2f
    ldr     q24, [x1]
    ldr     q25, [x2], #16
    ldr     q26, [x3]
    ldr     q27, [x4]
    fmul    v24.4s, v24.4s, v19.4s
    fmul    v26.4s, v26.4s, v16.4s
    fmla    v26.4s, v25.4s, v22.4s
    fmul    v27.4s, v27.4s, v17.4s
    fmul    v25.4s, v25.4s, v25.4s
    fmla    v27.4s, v25.4s, v23.4s
    fsqrt   v28.4s, v27.4s
    fmul    v28.4s, v28.4s, v21.4s
    fadd    v28.4s, v28.4s, v18.4s
    fmul    v29.4s, v26.4s, v20.4s
    fdiv    v29.4s, v29.4s, v28.4s
    fsub    v24.4s, v24.4s, v29.4s
    str     q24, [x1], #16
    str     q26, [x3], #16
    str     q27, [x4], #16
    sub     x0, x0, #4
    b       1b
2:  cbz     x0, 3f
    ldr     s24, [x1]
    ldr     s25, [x2], #4
    ldr     s26, [x3]
    ldr     s27, [x4]
    fmul    s24, s24, s19
    fmul    s26, s26, s16
    fmadd   s26, s25, s22, s26
    fmul    s27, s27, s17
    fmul    s25, s25, s25
    fmadd   s27, s25, s23, s27
    fsqrt   s28, s27
    fmul    s28, s28, s21
    fadd    s28, s28, s18
    fmul    s29, s26, s20
    fdiv    s29, s29, s28
    fsub    s24, s24, s29
    str     s24, [x1], #4
    str     s26, [x3], #4
    str     s27, [x4], #4
    sub     x0, x0, #1
    b       2b
3:  ret
