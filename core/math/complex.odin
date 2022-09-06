package math

// Trigonometric functions
cacos :: proc{cacos_d, cacos_f}
casin :: proc{casin_d, casin_f}
catan :: proc{catan_d, catan_f}
ccos  :: proc{ccos_d, ccos_f}
csin  :: proc{csin_d, csin_f}
ctan  :: proc{ctan_d, ctan_f}

cacos_d :: proc "contextless" (z: complex128) -> complex128 {
	z := casin(z)
	return complex(PI/2 - real(z), -imag(z))
}
cacos_f :: proc "contextless" (z: complex64) -> complex64 {
	z := casin(z)
	return complex(PI/2 - real(z), -imag(z))
}

casin_d :: proc "contextless" (z: complex128) -> complex128 {
	x := real(z)
	y := imag(z)
	w := complex(1.0 - (x - y)*(x + y), -2.0*x*y)
	r := clog(complex(-y, x) + csqrt(w))
	return complex(imag(r), -real(r))
}

casin_f :: proc "contextless" (z: complex64) -> complex64 {
	x := real(z)
	y := imag(z)
	w := complex(1.0 - (x - y)*(x + y), -2.0*x*y)
	r := clog(complex(-y, x) + csqrt(w))
	return complex(imag(r), -real(r))
}

/* origin: OpenBSD /usr/src/lib/libm/src/s_catan.c */
/*
 * Copyright (c) 2008 Stephen L. Moshier <steve@moshier.net>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */
/*
 *      Complex circular arc tangent
 *
 *
 * SYNOPSIS:
 *
 * double complex catan();
 * double complex z, w;
 *
 * w = catan (z);
 *
 *
 * DESCRIPTION:
 *
 * If
 *     z = x + iy,
 *
 * then
 *          1       (    2x     )
 * Re w  =  - arctan(-----------)  +  k PI
 *          2       (     2    2)
 *                  (1 - x  - y )
 *
 *               ( 2         2)
 *          1    (x  +  (y+1) )
 * Im w  =  - log(------------)
 *          4    ( 2         2)
 *               (x  +  (y-1) )
 *
 * Where k is an arbitrary integer.
 *
 * catan(z) = -i catanh(iz).
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       -10,+10      5900       1.3e-16     7.8e-18
 *    IEEE      -10,+10     30000       2.3e-15     8.5e-17
 * The check catan( ctan(z) )  =  z, with |x| and |y| < PI/2,
 * had peak relative error 1.5e-16, rms relative error
 * 2.9e-17.  See also clog().
 */

 // TODO: musl catan.c had this unused define. Important?:
 // #define MAXNUM 1.0e308
_redupi :: proc "contextless" (x: f64) -> f64 {
	DP1 :: f64(3.14159265160560607910E0)
	DP2 :: f64(1.98418714791870343106E-9)
	DP3 :: f64(1.14423774522196636802E-17)

	t := x/PI
	if t >= 0.0 {
		t += 0.5
	} else {
		t -= 0.5
	}

	// TODO: verify these two lines. implicit casts in the C impl
	i := i64(t)  /* the multiple */
	t = f64(i)
	t = ((x - t * DP1) - t * DP2) - t * DP3
	return t
}

catan_d :: proc "contextless" (z: complex128) -> complex128 {
	x := real(z)
	y := imag(z)

	x2 := x * x
	a := 1.0 - x2 - (y * y)

	t := 0.5 * atan2(2.0 * x, a)
	// TODO: verify that this line and the line before the return are equivalent. musl catan.h has an implicit cast here to
	// complex128
	wr := _redupi(t)

	t = y - 1.0
	a = x2 + (t * t)

	t = y + 1.0
	a = (x2 + t * t)/a
	w := complex(wr, 0.25 * ln(a))
	return w
}

/* origin: OpenBSD /usr/src/lib/libm/src/s_catanf.c */
/*
 * Copyright (c) 2008 Stephen L. Moshier <steve@moshier.net>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */
/*
 *      Complex circular arc tangent
 *
 *
 * SYNOPSIS:
 *
 * float complex catanf();
 * float complex z, w;
 *
 * w = catanf( z );
 *
 *
 * DESCRIPTION:
 *
 * If
 *     z = x + iy,
 *
 * then
 *          1       (    2x     )
 * Re w  =  - arctan(-----------)  +  k PI
 *          2       (     2    2)
 *                  (1 - x  - y )
 *
 *               ( 2         2)
 *          1    (x  +  (y+1) )
 * Im w  =  - log(------------)
 *          4    ( 2         2)
 *               (x  +  (y-1) )
 *
 * Where k is an arbitrary integer.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -10,+10     30000        2.3e-6      5.2e-8
 */

// TODO: musl catan.c had this unused define. Important?:
// #define MAXNUMF 1.0e38F

_redupif :: proc "contextless" (xx: f32) -> f32 {
	DP1 :: 3.140625
	DP2 :: 9.67502593994140625E-4
	DP3 :: 1.509957990978376432E-7

	x := xx
	t := x/PI
	if t >= 0.0 {
		t += 0.5
	} else {
		t -= 0.5
	}

	i := i64(t)  /* the multiple */
	t = f32(i)
	t = ((x - t * DP1) - t * DP2) - t * DP3
	return t
}

catan_f :: proc "contextless" (z: complex64) -> complex64 {
	x := real(z)
	y := imag(z)

	x2 := x * x
	a := 1.0 - x2 - (y * y)

	t := 0.5 * atan2(2.0 * x, a)
	wr := _redupif(t)

	t = y - 1.0
	a = x2 + (t * t)

	t = y + 1.0
	a = (x2 + (t * t))/a
	w := complex(wr, 0.25 * ln(a))
	return w
}

ccos_d :: proc "contextless" (z: complex128) -> complex128 {
	return ccosh(complex(-imag(z), real(z)))
}
ccos_f :: proc "contextless" (z: complex64) -> complex64 {
	return ccosh(complex(-imag(z), real(z)))
}
csin_d :: proc "contextless" (z: complex128) -> complex128 {
	z := csinh(complex(-imag(z), real(z)))
	return complex(imag(z), -real(z))
}
csin_f :: proc "contextless" (z: complex64) -> complex64 {
	z := csinh(complex(-imag(z), real(z)))
	return complex(imag(z), -real(z))
}
ctan_d :: proc "contextless" (z: complex128) -> complex128 {
	z := ctanh(complex(-imag(z), real(z)))
	return complex(imag(z), -real(z))
}
ctan_f :: proc "contextless" (z: complex64) -> complex64 {
	z := ctanh(complex(-imag(z), real(z)))
	return complex(imag(z), -real(z))
}

// Hyperbolic functions
cacosh :: proc{cacosh_d, cacosh_f}
casinh :: proc{casinh_d, casinh_f}
catanh :: proc{catanh_d, catanh_f}
ccosh  :: proc{ccosh_d, ccosh_f}
csinh  :: proc{csinh_d, csinh_f}
ctanh  :: proc{ctanh_d, ctanh_f}

cacosh_d :: proc "contextless" (z: complex128) -> complex128 {
	zineg := signbit(imag(z))

	z := cacos(z)
	if zineg {
		return complex(imag(z), -real(z))
	} else {
		return complex(-imag(z), real(z))
	}
}
cacosh_f :: proc "contextless" (z: complex64) -> complex64 {
	zineg := signbit(imag(z))

	z := cacos(z)
	if zineg {
		return complex(imag(z), -real(z))
	} else {
		return complex(-imag(z), real(z))
	}
}
casinh_d :: proc "contextless" (z: complex128) -> complex128 {
	z := casin(complex(-imag(z), real(z)))
	return complex(imag(z), -real(z))
}
casinh_f :: proc "contextless" (z: complex64) -> complex64 {
	z := casin(complex(-imag(z), real(z)))
	return complex(imag(z), -real(z))
}
catanh_d :: proc "contextless" (z: complex128) -> complex128 {
	z := catan(complex(-imag(z), real(z)))
	return complex(imag(z), -real(z))
}
catanh_f :: proc "contextless" (z: complex64) -> complex64 {
	z := catan(complex(-imag(z), real(z)))
	return complex(imag(z), -real(z))
}

ccosh_d :: proc "contextless" (z: complex128) -> complex128 {
	// TODO: complicated
	return {}
}
ccosh_f :: proc "contextless" (z: complex64) -> complex64 {
	// TODO: complicated
	return {}
}
csinh_d :: proc "contextless" (z: complex128) -> complex128 {
	// TODO: complicated
	return {}
}
csinh_f :: proc "contextless" (z: complex64) -> complex64 {
	// TODO: complicated
	return {}
}
ctanh_d :: proc "contextless" (z: complex128) -> complex128 {
	// TODO: complicated
	return {}
}
ctanh_f :: proc "contextless" (z: complex64) -> complex64 {
	// TODO: complicated
	return {}
}

// Exponential and logarithmic functions
cexp :: proc{cexp_d, cexp_f}
clog :: proc{clog_d, clog_f}

cexp_d :: proc "contextless" (z: complex128) -> complex128 {
	// TODO: complicated
	return {}
}

cexp_f :: proc "contextless" (z: complex64) -> complex64 {
	// TODO: complicated
	return {}
}

clog_d :: proc "contextless" (z: complex128) -> complex128 {
	return {}
}

clog_f :: proc "contextless" (z: complex64) -> complex64 {
	return {}
}

// Power and absolute-value functions
cabs  :: proc{cabs_d, cabs_f}
cpow  :: proc{cpow_d, cpow_f}
csqrt :: proc{csqrt_d, csqrt_f}

cabs_d :: proc "contextless" (z: complex128) -> complex128 {
	return {}
}

cabs_f :: proc "contextless" (z: complex64) -> complex64 {
	return {}
}

cpow_d :: proc "contextless" (z: complex128) -> complex128 {
	return {}
}

cpow_f :: proc "contextless" (z: complex64) -> complex64 {
	return {}
}

csqrt_d :: proc "contextless" (z: complex128) -> complex128 {
	return {}
}

csqrt_f :: proc "contextless" (z: complex64) -> complex64 {
	return {}
}

// Manipulation functions
carg  :: proc{carg_d, carg_f}
/*cimag :: proc{cimag_d, cimag_f}*/
conj  :: proc{conj_d, conj_f}
cproj :: proc{cproj_d, cproj_f}
/*creal :: proc{creal_d, creal_f}*/

carg_d :: proc "contextless" (z: complex128) -> f64 {
	return {}
}

carg_f :: proc "contextless" (z: complex64) -> f32 {
	return {}
}

/*cimag_d :: proc "contextless" (z: complex128) -> f64 {*/
	/*return {}*/
/*}*/

/*cimag_f :: proc "contextless" (z: complex64) -> f32 {*/
	/*return {}*/
/*}*/

conj_d :: proc "contextless" (z: complex128) -> complex128 {
	return {}
}

conj_f :: proc "contextless" (z: complex64) -> complex64 {
	return {}
}

cproj_d :: proc "contextless" (z: complex128) -> complex128 {
	return {}
}

cproj_f :: proc "contextless" (z: complex64) -> complex64 {
	return {}
}

/*creal_d :: proc "contextless" (z: complex128) -> f64 {*/
	/*return {}*/
/*}*/

/*creal_f :: proc "contextless" (z: complex64) -> f32 {*/
	/*return {}*/
/*}*/
