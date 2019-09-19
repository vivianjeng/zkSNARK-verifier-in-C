
#include <stdio.h>
#include <gmp.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#ifndef CPUCYCLES_amd64cpuinfo_h
#define CPUCYCLES_amd64cpuinfo_h

#ifdef __cplusplus
extern "C" {
#endif

extern long long cpucycles_amd64cpuinfo(void);
extern long long cpucycles_amd64cpuinfo_persecond(void);

#ifdef __cplusplus
}
#endif

#ifndef cpucycles_implementation
#define cpucycles_implementation "amd64cpuinfo"
#define cpucycles cpucycles_amd64cpuinfo
#define cpucycles_persecond cpucycles_amd64cpuinfo_persecond
#endif

#endif

#if(GMP_LIMB_BITS == 32)
#define N_LIMBS 8
#elif(GMP_LIMB_BITS == 64)
#define N_LIMBS 4
#else
#error "Only 32 and 64 bit architectures are supported"
#endif

// Parameters used in fp.c
#define BN_P "21888242871839275222246405745257275088696311157297823662689037894645226208583"
#define BN_PINV32 3834012553UL
#define BN_PINV64 9786893198990664585UL

// Parameters used in fp2.c
#define ALPHA (-1) // constant coefficient in the irreducible polynomial x^2 - alpha, used to construct F_{p^2}

// Parameters used in curve.c
#define BN_X "4965661367192848881" // parameter x used to generate the curve (see "Pairing-Friendly Elliptic Curves of Prime Order")
#define BN_N "21888242871839275222246405745257275088548364400416034343698204186575808495617" // prime order of E(F_p)
#define BN_TRACE "147946756881789318990833708069417712967" // trace of Frobenius of the curve
#define BN_CHI "552484233613224096312617126783173147097382103762957654188882734314196910839907541213974502761540629817009608548654680343627701153829446747810907373256841551006201639677726139946029199968412598804882391702273019083653272047566316584365559776493027495458238373902875937659943504873220554161550525926302303331747463515644711876653177129578303191095900909191624817826566688241804408081892785725967931714097716709526092261278071952560171111444072049229123565057483750161460024353346284167282452756217662335528813519139808291170539072125381230815729071544861602750936964829313608137325426383735122175229541155376346436093930287402089517426973178917569713384748081827255472576937471496195752727188261435633271238710131736096299798168852925540549342330775279877006784354801422249722573783561685179618816480037695005515426162362431072245638324744480" // exponent of final exponentiation
#define BN_CHICOMP "4965661367192848881" // exponent of final exponentiation for compressed pairing
#define BN_LOOPLENGTH_ETA "11916685325773193009570696613837024235910666865622157670373"
#define BN_B "3" // parameter b in the curve equation x^2 = y^3 + b

// Parameters used in fp6.c
#define BN_XI "1", "9"
#define BN_YPMINUS1 "10307601595873709700152284273816112264069230130616436755625194854815875713954", "21575463638280843010398324269430826099269044274347216827212613867836435027261"
#define BN_ZETA "2203960485148121921418603742825762020974279258880205651966"
// #define BN_XI2 "2", "0"
// #define BN_1O27XI3 "5674729633439812094656475563585219467439784374114250579215676491204317905929", "16213513238399463127589930181672055621256526783183573083473361403440908302654"
// #define BN_1O3XI3 "7296080957279758407415468581752425029565437052432607887563012631548408736195", "14592161914559516814830937163504850059130874104865215775126025263096817472388"
// #define BN_1O3XI "14592161914559516814830937163504850059130874104865215775126025263096817472389", "14592161914559516814830937163504850059130874104865215775126025263096817472389"
// #define BN_1O3MODP "14592161914559516814830937163504850059130874104865215775126025263096817472389"
// #define BN_COMETA_C0_CONST "82434016654300679717245125061265641166427769693017678351449950981238451124296"
// #define BN_COMETA_C1_CONST "82434016654300679719231239282227840001499775752201953636308636697028740308739"

// Parameters used in fp12.c
#define BN_TAU "0", "0", "0", "1", "0", "0" // constant tau used to construct F_p^12 as F_p^6[Z]/ (Z^2 - tau)
#define BN_ZPMINUS1 "16469823323077808223889137241176536799009286646108169935659301613961712198316", "8376118865763821496583973867626364092589906065868298776909617916018768340080" // Z^(p-1)
#define BN_ZPMINUS1INV "5722266937896532885780051958958348231143373700109372999374820235121374419868", "18566938241244942414004596690298913868373833782006617400804628704885040364344" // Z^(1-p)

// Parameters used in points.c
#define BN_CURVEGEN "1", "2", "1"
#define BN_TWISTGEN_X "11559732032986387107991004021392285783925812861821192530917403151452391805634", "10857046999023057135944570762232829481370756359578518086990519993285655852781"
#define BN_TWISTGEN_Y "4082367875863433681332203403145435568316851327593401208105741076214120093531", "8495653923123431417604973247489272438418190587263600148770280649306958101930"

// Parameters used for OptAte computation in ate_optate.c
#define BN_ZETA2 "21888242871839275220042445260109153167277707414472061641714758635765020556616" // zeta^2
#define BN_Z2P "10307601595873709700152284273816112264069230130616436755625194854815875713954", "21575463638280843010398324269430826099269044274347216827212613867836435027261" // Z^(2p)
#define BN_Z3P "3505843767911556378687030309984248845540243509899259641013678093033130930403", "2821565182194536844548159561693502659359617185244120367078079554186484126554" // Z^(3p)

mpz_t p;
unsigned long p_inv; // -p^{-1} mod 2^{GMP_LIMB_BITS} used in Montgomery reduction

typedef struct fpe_struct fpe_struct_t;

struct fpe_struct
{
	mp_limb_t m_value[N_LIMBS];
};

void fp_init()
{
	mpz_init_set_str(p, BN_P, 10);

	// Set the value m', used in Montgomery reduction, see HAC, alg. 14.32
#if(GMP_LIMB_BITS == 32)
	p_inv = BN_PINV32;
#elif(GMP_LIMB_BITS == 64)
	p_inv = BN_PINV64;
#else
#error "Only 32 and 64 Bit architectures are supported"
#endif
}

typedef fpe_struct_t fpe_t[1];

// Montgomery reduction, see HAC, alg 14.32, size of op must be 2 * N_LIMBS:
static void fpe_montgomery_reduce(fpe_t rop, mp_limb_t *op)
{
	mp_limb_t u, t, c, d; // dummies
	int i;
	c = 0;

	for(i = 0; i < N_LIMBS; i++)
	{
		u = op[i] * p_inv;
		d = mpn_addmul_1(op + i, p->_mp_d, N_LIMBS, u);
		t = op[N_LIMBS + i] + d;
		op[N_LIMBS + i] = t + c;
		if(t < d || (c == 1 && t + c == 0))
			c = 1;
		else
			c = 0;
	}
	if(c || mpn_cmp(op + N_LIMBS, p->_mp_d, N_LIMBS) >= 0)
		// Result is larger than p, subtract p from the result:
		mpn_sub_n(rop->m_value, op + N_LIMBS, p->_mp_d, N_LIMBS);
	else
	{
		for(i = 0; i < N_LIMBS; i++)
			rop->m_value[i] = op[N_LIMBS + i];
	}
}

// Montgomery transformation:
static void fpe_montgomery_trans(fpe_t rop)
{
	mp_limb_t tmp1[2 * N_LIMBS], tmp2[2 * N_LIMBS]; // needed for intermediate results
	int i;

	for(i = 0; i < N_LIMBS; i++)
		tmp1[i] = 0;
	for(i = 0; i < N_LIMBS; i++)
		tmp1[N_LIMBS + i] = rop->m_value[i];
	mpn_tdiv_qr(tmp2, rop->m_value, 0, tmp1, 2 * (N_LIMBS), p->_mp_d, N_LIMBS);
}

// Montgomery retransformation:
static void fpe_montgomery_retrans(fpe_t rop, const fpe_t op)
{
	mp_limb_t tmp[2 * N_LIMBS]; // needed for intermediate results
	int i;

	for(i = 0; i < N_LIMBS; i++)
		tmp[N_LIMBS + i] = 0;
	for(i = 0; i < N_LIMBS; i++)
		tmp[i] = op->m_value[i];
	fpe_montgomery_reduce(rop, tmp);
}

void fpe_set(fpe_t rop, const fpe_t op)
{
	int i;
	for(i = 0; i < N_LIMBS; i++)
		rop->m_value[i] = op->m_value[i];
}

// Set fpe_t rop to given value:
void fpe_set_ui(fpe_t rop, const unsigned int op)
{
	int i;
	for(i = 1; i < N_LIMBS; i++)
		rop->m_value[i] = 0;
	rop->m_value[0] = op;

	fpe_montgomery_trans(rop);
}

// Set fpe_t rop to given value given as (ASCII) string
void fpe_set_str(fpe_t rop, const char* op)
{

	// Initialize all limbs with 0:
	int i;
	for(i = 0; i < N_LIMBS; i++)
	{
		rop->m_value[i] = 0;
	}
	// Determine the length of op:
	const char *scan = op;
	int size = 0;
	while(*scan != 0)
	{
		++scan;
		++size;
	}

	unsigned char str[size];
	memcpy(str, op, size);

	// Convert from ASCII:
	for(i = 0; i < size; ++i)
		str[i] -= '0';

	mpn_set_str(rop->m_value, str, size, 10);

	fpe_montgomery_trans(rop);
}

// Set rop to one
void fpe_setone(fpe_t rop)
{
	int i;
	for(i = 1; i < N_LIMBS; i++)
		rop->m_value[i] = 0;
	rop->m_value[0] = 1;
	fpe_montgomery_trans(rop);
}

// Set rop to zero
void fpe_setzero(fpe_t rop)
{
	int i;
	for(i = 0; i < N_LIMBS; i++)
		rop->m_value[i] = 0;
}

// Return 1 if op is zero, 0 otherwise
int fpe_iszero(const fpe_t op)
{
	int i;
	for(i = 0; i < N_LIMBS; i++)
		if(op->m_value[i]) return 0;
	return 1;
}

// Return 1 if op is one, 0 otherwise
int fpe_isone(const fpe_t op)
{
	int i;

	for(i = 1; i < N_LIMBS; i++)
		if(!op->m_value[i]) return 0;

	// Retransform from Montgomery domain:
	fpe_t dummy;
	fpe_montgomery_retrans(dummy, op);

	return dummy->m_value[0] == 1;
}

// Compute the negative of an fpe
void fpe_neg(fpe_t rop, const fpe_t op)
{
	mpn_sub_n(rop->m_value, p->_mp_d, op->m_value, N_LIMBS);
}

// Double an fpe:
void fpe_double(fpe_t rop, const fpe_t op)
{
	mp_limb_t c;
	c = mpn_lshift(rop->m_value, op->m_value, N_LIMBS, 1);
	
	// Reduce if result is larger than p: 
	if(c || mpn_cmp(rop->m_value, p->_mp_d, N_LIMBS) > 0)
		mpn_sub_n(rop->m_value, rop->m_value, p->_mp_d, N_LIMBS);
}

// Halve an fpe:
void fpe_halve(fpe_t rop, const fpe_t op)
{
	if((op->m_value[0] % 2) == 0)
		mpn_rshift(rop->m_value, op->m_value, N_LIMBS, 1);
	else
	{
		int c;
		c = (mpn_add_n(rop->m_value, op->m_value, p->_mp_d, N_LIMBS)) << (GMP_LIMB_BITS - 1);
		mpn_rshift(rop->m_value, rop->m_value, N_LIMBS, 1);
		rop->m_value[N_LIMBS - 1] ^= c;
	}
}

// Triple an fpe:
void fpe_triple(fpe_t rop, const fpe_t op)
{
	mp_limb_t tmp[2 * N_LIMBS]; // needed for intermediate results
	int i, c;

	for(i = 2 * N_LIMBS - 1; i >= N_LIMBS; i--)
		tmp[i] = 0;

	// Double
	c = mpn_lshift(tmp, op->m_value, N_LIMBS, 1);
	if(c || mpn_cmp(tmp, p->_mp_d, N_LIMBS) >= 0)
		mpn_sub_n(tmp, tmp, p->_mp_d, N_LIMBS);

	// Add
	c = mpn_add_n(rop->m_value, tmp, op->m_value, N_LIMBS);	
	if(c || mpn_cmp(rop->m_value, p->_mp_d, N_LIMBS) >= 0)
		mpn_sub_n(rop->m_value, rop->m_value, p->_mp_d, N_LIMBS);
}

// Add two fpe, store result in rop:
void fpe_add(fpe_t rop, const fpe_t op1, const fpe_t op2)
{
	mp_limb_t c;
	c = mpn_add_n(rop->m_value, op1->m_value, op2->m_value, N_LIMBS);

	// Reduce if result is larger than p: 
	if(c || mpn_cmp(rop->m_value, p->_mp_d, N_LIMBS) >= 0)
		mpn_sub_n(rop->m_value, rop->m_value, p->_mp_d, N_LIMBS);
}

// Subtract op2 from op1, store result in rop:
void fpe_sub(fpe_t rop, const fpe_t op1, const fpe_t op2)
{
	mp_limb_t b;
	b = mpn_sub_n(rop->m_value, op1->m_value, op2->m_value, N_LIMBS);
	
	if(b)
		mpn_add_n(rop->m_value, rop->m_value, p->_mp_d, N_LIMBS);
}

// Multiply two fpe, store result in rop:
void fpe_mul(fpe_t rop, const fpe_t op1, const fpe_t op2)
{
#ifdef BENCH
    nummultp++;
    multpcycles -= cpucycles();
#endif

	mp_limb_t tmp[2 * N_LIMBS]; // needed for intermediate results


	if(fpe_iszero(op1) || fpe_iszero(op2))
		fpe_setzero(rop);
	else
	{
		mpn_mul_n(tmp, op1->m_value, op2->m_value, N_LIMBS);
		fpe_montgomery_reduce(rop, tmp);
	}

#ifdef BENCH
    multpcycles += cpucycles();
#endif
}

// Compute inverse of an fpe, store result in rop:
void fpe_invert(fpe_t rop, const fpe_t op)
{
#ifdef BENCH
    numinvp++;
    invpcycles -= cpucycles();
#endif
	/*
	 * FIXME: This code doesn't work they way it should
	// Using mpn_gcdext:
	mp_limb_t mp_limb_t d1[N_LIMBS + 1], mp_limb_t d2[N_LIMBS + 1], tmp1[2 * N_LIMBS], tmp2[2 * N_LIMBS]; // needed for intermediate results
	mp_size_t n1, n2;
	
	int i;
	for(i = 0; i < N_LIMBS; i++)
	{
		d1[i] = p->_mp_d[i];
		d2[i] = op->m_value[i];
	}

	n2 = mpn_gcdext(tmp1, tmp2, &n1, d2, N_LIMBS, d1, N_LIMBS);
	

	for(i = 0; i < N_LIMBS; i++)
		rop->m_value[i] = tmp2[i];

	fpe_montgomery_trans(rop);
	fpe_montgomery_trans(rop);
	
	*/


	// Partial Montgomery inversion, see Guide to ECC, alg. 2.23
	mp_limb_t u[N_LIMBS], v[N_LIMBS], s[N_LIMBS];
	unsigned long i, k, c;
	c = k = 0; 
	for(i = 0; i < N_LIMBS; i++)
	{
		u[i] = p->_mp_d[i];
		v[i] = op->m_value[i];
		rop->m_value[i] = 0;
		s[i] = 0;
	}
	s[0] = 1;

	int v_zero = 1;
	for(i = 0; i < N_LIMBS && v_zero; i++)
		if(v[i]) v_zero = 0;
	// Don't try to invert 0:
	assert(!v_zero);

	while(!v_zero)
	{
		if(u[0] % 2 == 0)
		{
			mpn_rshift(u, u, N_LIMBS, 1);
			mpn_lshift(s, s, N_LIMBS, 1);
		}
		else
		{
			if(v[0] % 2 == 0)
			{
				mpn_rshift(v, v, N_LIMBS, 1);
				mpn_lshift(rop->m_value, rop->m_value, N_LIMBS, 1);
			}
			else
			{
				if(mpn_cmp(u, v, N_LIMBS) > 0)
				{
					mpn_sub_n(u, u, v, N_LIMBS);
					mpn_rshift(u, u, N_LIMBS, 1);
					c = mpn_add_n(rop->m_value, rop->m_value, s, N_LIMBS);
					mpn_lshift(s, s, N_LIMBS, 1);
				}
				else
				{
					mpn_sub_n(v, v, u, N_LIMBS);
					mpn_rshift(v, v, N_LIMBS, 1);
					mpn_add_n(s, s, rop->m_value, N_LIMBS);
					c = mpn_lshift(rop->m_value, rop->m_value, N_LIMBS, 1);
				}
			}
		}
		++k;
		// Check, whether v is zero:
		v_zero = 1;
		for(i = 0; i < N_LIMBS; i++)
		{
			if(v[i]) v_zero = 0;
		}
	}
	
	if(c || mpn_cmp(rop->m_value, p->_mp_d, N_LIMBS) >= 0)
		mpn_sub_n(rop->m_value, rop->m_value, p->_mp_d, N_LIMBS);

	mpn_sub_n(rop->m_value, p->_mp_d, rop->m_value, N_LIMBS);
	
	// Make the Montgomery Inversion complete:
	for(; k < N_LIMBS * sizeof(mp_limb_t) * 16; k++)
	{
		c = mpn_lshift(rop->m_value, rop->m_value, N_LIMBS, 1);
		if(c || mpn_cmp(rop->m_value, p->_mp_d, N_LIMBS) >= 0)
			mpn_sub_n(rop->m_value, rop->m_value, p->_mp_d, N_LIMBS);
	}
#ifdef BENCH
    invpcycles += cpucycles();
#endif
}

// Print the element to stdout:
void fpe_print(FILE *outfile, const fpe_t op)
{
	int i;

	// Retransform from Montgomery domain:
	fpe_t dummy;
	fpe_montgomery_retrans(dummy, op);

	i = 0;
	
	while(i < N_LIMBS && !dummy->m_value[N_LIMBS - 1 - i])
		i++;
	// Print '0' if op is zero:
	if(i == N_LIMBS && !dummy->m_value[0])
	{
		fputc('0', outfile);
		fputc(0, outfile);
	}
	else
	{
		unsigned char *str;
		double log2 = M_LOG10E * M_LN2; // log_{10}2
		
		size_t str_size = log2 * (GMP_LIMB_BITS * (N_LIMBS)) + 2;
		str = malloc(str_size);

		str_size = mpn_get_str(str, 10, dummy->m_value, N_LIMBS - i);
	
		// Strip leading zeros:
		while(*str == 0 && str_size > 1)
		{
			str++;
			str_size--;
		}
		
		for(i = 0; i < str_size; i++)
		{
			fputc(*(str + i) + '0', outfile);
		}
		fputc(0, outfile);
		free(str);
	}
}

/// Structure describing a point on a BN-curve
typedef struct curvepoint_fp_struct curvepoint_fp_struct_t;
struct curvepoint_fp_struct
{	
	fpe_t m_x; // X-Coordinate (Jacobian Coordinate system)
	fpe_t m_y; // Y-Coordinate (Jacobian Coordinate system)
	fpe_t m_z; // Y-Coordinate (Jacobian Coordinate system)
	fpe_t m_t; // T = Z^2, only used during pairing computation, set to zero if not set
};

typedef curvepoint_fp_struct_t curvepoint_fp_t[1];

#define fpe_square(rop, op) fpe_mul(rop, op, op)

mpz_t b; /* parameter b in the curve equation y^2 = x^3 + b */
mpz_t n; /* order of the curve */
mpz_t x;  
mpz_t trace; /* trace of Frobenius */
mpz_t chi; /* p^12 / n */
mpz_t chicomp; 
mpz_t looplength_eta;
mpz_t ate_loop_count; 

void curve_init()
{
	/* Curve parameters */
	mpz_init_set_str(x, BN_X, 10);
	mpz_init_set_str(n, BN_N, 10);
	mpz_init_set_str(trace, BN_TRACE,10);
	mpz_init_set_str(chi, BN_CHI, 10); // (p^k - 1) / n
	mpz_init_set_str(chicomp, BN_CHICOMP, 10);
	mpz_init_set_str(looplength_eta, BN_LOOPLENGTH_ETA, 10);
	mpz_init_set_str(b, BN_B, 10);
}

// Global dummies usable by all curvepoints:
fpe_t curvepoint_dummy_fpe1;


void curvepoint_fp_init_set_str(curvepoint_fp_t rop, const char* x, const char* y, const char* z)
{
	fpe_set_str(rop->m_x, x);
	fpe_set_str(rop->m_y, y);
	fpe_set_str(rop->m_z, z);
	fpe_set_ui(rop->m_t, 0);
}

// Set the coordinates of a curvepoint_fp_t by copying the coordinates from another curvepoint_fp
void curvepoint_fp_set(curvepoint_fp_t rop, const curvepoint_fp_t op)
{
	fpe_set(rop->m_x, op->m_x);
	fpe_set(rop->m_y, op->m_y);
	fpe_set(rop->m_z, op->m_z);
	fpe_set_ui(rop->m_t, 0);
}

// Set the coordinates of a curvepoint_fp:
void curvepoint_fp_set_str(curvepoint_fp_t rop, const char* x, const char* y, const char* z)
{
	fpe_set_str(rop->m_x, x);
	fpe_set_str(rop->m_y, y);
	fpe_set_str(rop->m_z, z);
	fpe_set_ui(rop->m_t, 0);
}

// Addition of two points, op2 is assumed to be in affine coordinates 
// For the algorithm see e.g. DA Peter Schwabe
void curvepoint_fp_mixadd(curvepoint_fp_t rop, const curvepoint_fp_t op1, const curvepoint_fp_t op2)
{
	fpe_t tfpe1, tfpe2, tfpe3, tfpe4, tfpe5, tfpe6, tfpe7, tfpe8, tfpe9; // Temporary variables needed for intermediary results
	fpe_square(tfpe1, op1->m_z);
	fpe_mul(tfpe2, op1->m_z, tfpe1);
	fpe_mul(tfpe3, op2->m_x, tfpe1);
	fpe_mul(tfpe4, op2->m_y, tfpe2);
	fpe_sub(tfpe5, tfpe3, op1->m_x);
	fpe_sub(tfpe6, tfpe4, op1->m_y);
	fpe_square(tfpe7, tfpe5);
	fpe_mul(tfpe8, tfpe7, tfpe5);
	fpe_mul(tfpe9, op1->m_x, tfpe7);

	fpe_double(tfpe1, tfpe9);
	fpe_add(tfpe1, tfpe1, tfpe8);
	fpe_square(rop->m_x, tfpe6);
	fpe_sub(rop->m_x, rop->m_x, tfpe1);
	fpe_sub(tfpe1, tfpe9, rop->m_x);
	fpe_mul(tfpe2, tfpe1, tfpe6);
	fpe_mul(tfpe3, op1->m_y, tfpe8);
	fpe_sub(rop->m_y, tfpe2, tfpe3);
	fpe_mul(rop->m_z, op1->m_z, tfpe5);
}

void curvepoint_fp_double(curvepoint_fp_t rop, const curvepoint_fp_t op)
{
	fpe_t tfpe1, tfpe2, tfpe3, tfpe4; // Temporary variables needed for intermediary results
	fpe_square(tfpe1, op->m_y);
	fpe_mul(tfpe2, tfpe1, op->m_x);
	fpe_double(tfpe2, tfpe2);
	fpe_double(tfpe2, tfpe2);
	fpe_square(tfpe3, tfpe1);
	fpe_double(tfpe3, tfpe3);
	fpe_double(tfpe3, tfpe3);
	fpe_double(tfpe3, tfpe3);
	fpe_square(tfpe4, op->m_x);
	fpe_triple(tfpe4, tfpe4);
	fpe_square(rop->m_x, tfpe4);
	fpe_double(tfpe1, tfpe2);
	fpe_sub(rop->m_x, rop->m_x, tfpe1);
	fpe_sub(tfpe1, tfpe2, rop->m_x);
	fpe_mul(rop->m_z, op->m_y, op->m_z);
	fpe_double(rop->m_z, rop->m_z);
	fpe_mul(rop->m_y, tfpe4, tfpe1);
	fpe_sub(rop->m_y, rop->m_y, tfpe3);
}

void curvepoint_fp_mul(curvepoint_fp_t rop, const curvepoint_fp_t op, const mpz_t scalar)
{
	size_t i;
	curvepoint_fp_t r;
	curvepoint_fp_set(r, op);
	for(i = mpz_sizeinbase(scalar, 2) - 1; i > 0; i--)
	{
		curvepoint_fp_double(r, r);
		if(mpz_tstbit(scalar, i - 1)) 
			curvepoint_fp_mixadd(r, r, op);
	}
	curvepoint_fp_set(rop, r);
}



// Negate a point, store in rop:
void curvepoint_fp_neg(curvepoint_fp_t rop, const curvepoint_fp_t op)
{
	fpe_neg(curvepoint_dummy_fpe1, op->m_y);
	fpe_set(rop->m_x, op->m_x);
	fpe_set(rop->m_y, curvepoint_dummy_fpe1);
	fpe_set(rop->m_z, op->m_z);
}

// Transform to Affine Coordinates (z=1)
void curvepoint_fp_makeaffine(curvepoint_fp_t point)
{
	if(fpe_iszero(point->m_z))
	{
		fpe_setzero(point->m_x);
		fpe_setone(point->m_y);
		fpe_setzero(point->m_z);
	}
	else
	{
		fpe_invert(curvepoint_dummy_fpe1, point->m_z);
		fpe_mul(point->m_x, point->m_x, curvepoint_dummy_fpe1);
		fpe_mul(point->m_x, point->m_x, curvepoint_dummy_fpe1);

		fpe_mul(point->m_y, point->m_y, curvepoint_dummy_fpe1);
		fpe_mul(point->m_y, point->m_y, curvepoint_dummy_fpe1);
		fpe_mul(point->m_y, point->m_y, curvepoint_dummy_fpe1);

		fpe_setone(point->m_z);
	}
}

// Print a point:
void curvepoint_fp_print(FILE *outfile, const curvepoint_fp_t point)
{
	fprintf(outfile, "[");
	fpe_print(outfile, point->m_x);
	fprintf(outfile, ", ");
	fpe_print(outfile, point->m_y);
	fprintf(outfile, ", ");
	fpe_print(outfile, point->m_z);
	fprintf(outfile, "]");
}

// Elements from F_{p^2}= F_p[X] / (x^2 - alpha)F_p[X] are represented as aX + b
typedef struct fp2e_struct fp2e_struct_t;

struct fp2e_struct
{
	fpe_t m_a;
	fpe_t m_b;
};

typedef fp2e_struct_t fp2e_t[1];

// Set fp2e_t rop to given value:
void fp2e_set(fp2e_t rop, const fp2e_t op)
{
	fpe_set(rop->m_a, op->m_a);
	fpe_set(rop->m_b, op->m_b);
}

// Set fp2e_t rop to given value:
void fp2e_set_fpe(fp2e_t rop, const fpe_t op)
{
	fpe_set_ui(rop->m_a, 0);
	fpe_set(rop->m_b, op);
}

// Set rop to one
void fp2e_setone(fp2e_t rop)
{
	fpe_setzero(rop->m_a);
	fpe_setone(rop->m_b);
}

// Set rop to zero
void fp2e_setzero(fp2e_t rop)
{
	fpe_setzero(rop->m_a);
	fpe_setzero(rop->m_b);;
}

// Set an fp2e_t to value given in two strings
void fp2e_set_str(fp2e_t rop, const char* a_str, const char* b_str)
{
	fpe_set_str(rop->m_a, a_str);
	fpe_set_str(rop->m_b, b_str);
}

// Double an fp2e:
void fp2e_double(fp2e_t rop, const fp2e_t op)
{
	fpe_double(rop->m_a, op->m_a);
	fpe_double(rop->m_b, op->m_b);
}

// Triple an fp2e:
void fp2e_triple(fp2e_t rop, const fp2e_t op)
{
	fpe_triple(rop->m_a, op->m_a);
	fpe_triple(rop->m_b, op->m_b);
}

// Add two fp2e, store result in rop:
void fp2e_add(fp2e_t rop, const fp2e_t op1, const fp2e_t op2)
{
	fpe_add(rop->m_a, op1->m_a, op2->m_a);
	fpe_add(rop->m_b, op1->m_b, op2->m_b);
}

// Subtract op2 from op1, store result in rop:
void fp2e_sub(fp2e_t rop, const fp2e_t op1, const fp2e_t op2)
{
	fpe_sub(rop->m_a, op1->m_a, op2->m_a);
	fpe_sub(rop->m_b, op1->m_b, op2->m_b);
}

// Negate op
void fp2e_neg(fp2e_t rop, const fp2e_t op)
{
	fpe_neg(rop->m_a, op->m_a);
	fpe_neg(rop->m_b, op->m_b);
}


// Multiply two fp2e, store result in rop:
void fp2e_mul(fp2e_t rop, const fp2e_t op1, const fp2e_t op2)
{
#ifdef BENCH
    nummultp2++;
    multp2cycles -= cpucycles();
#endif
	fpe_t tmp1, tmp2, tmp3; // Needed for intermediary results

	if((fpe_iszero(op1->m_a) && fpe_iszero(op1->m_b)) || (fpe_iszero(op2->m_a) && fpe_iszero(op2->m_b)))
		fp2e_setzero(rop);
	else
	{
		fpe_mul(tmp1, op1->m_a, op2->m_a);
		fpe_mul(tmp2, op1->m_b, op2->m_b);
		fpe_add(tmp3, op2->m_a, op2->m_b);
		fpe_add(rop->m_a, op1->m_a, op1->m_b);
		fpe_set(rop->m_b, tmp2);
		fpe_mul(rop->m_a, rop->m_a, tmp3);
		fpe_sub(rop->m_a, rop->m_a, tmp1);
		fpe_sub(rop->m_a, rop->m_a, rop->m_b);
#if(ALPHA == 2)
		fpe_double(tmp1, tmp1);
		fpe_add(rop->m_b, rop->m_b, tmp1);
#elif(ALPHA == -2)
		fpe_double(tmp1, tmp1);
		fpe_sub(rop->m_b, rop->m_b, tmp1);
#elif(ALPHA == -1)
		fpe_sub(rop->m_b, rop->m_b, tmp1);
#else
#error "ALPHA must be -1, 2 or -2"
#endif
	}
#ifdef BENCH
    multp2cycles += cpucycles();
#endif
}

// Square an fp2e, store result in rop:
void fp2e_square(fp2e_t rop, const fp2e_t op)
{
#ifdef BENCH
    numsqp2++;
    sqp2cycles -= cpucycles();
#endif
	fpe_t tmp1, tmp2, tmp3; // Needed for intermediary results

	fpe_mul(tmp1, op->m_a, op->m_b);

	fpe_add(tmp2, op->m_a, op->m_b);
#if(ALPHA == 2)
	fpe_double(tmp3, op->m_a);
	fpe_add(rop->m_b, op->m_b, tmp3);
#elif(ALPHA == -2)
	fpe_double(tmp3, op->m_a);
	fpe_sub(rop->m_b, op->m_b, tmp3);
#elif(ALPHA == -1)
	fpe_sub(rop->m_b, op->m_b, op->m_a);
#else
#error "ALPHA must be -1, 2 or -2"
#endif
	fpe_mul(rop->m_b, rop->m_b, tmp2);

	fpe_sub(rop->m_b, rop->m_b, tmp1);
#if(ALPHA == 2)
	fpe_double(tmp2, tmp1);
	fpe_sub(rop->m_b, rop->m_b, tmp2);
#elif(ALPHA == -2)
	fpe_double(tmp2, tmp1);
	fpe_add(rop->m_b, rop->m_b, tmp2);
#elif(ALPHA == -1)
	fpe_add(rop->m_b, rop->m_b, tmp1);
#else
#error "ALPHA must be -1, 2 or -2"
#endif

	fpe_double(rop->m_a, tmp1);
#ifdef BENCH
    sqp2cycles += cpucycles();
#endif
}

// Multiply by xi which is used to construct F_p^6
void fp2e_mulxi(fp2e_t rop, const fp2e_t op)
{
    //TODO Check for XI and ALPHA
    fpe_t tmp1, tmp2;
   
    fpe_sub(tmp1, op->m_b, op->m_a); 
    fpe_sub(tmp2, tmp1, op->m_a); 
    fpe_add(rop->m_a, op->m_b, op->m_a);
    fpe_set(rop->m_b, tmp2);

    // MODIFIED
    // fp2e_mul(rop, op, xi);
}

// Multiply an fpe by xi which is used to construct F_p^6
void fp2e_mulxi_fpe(fp2e_t rop, const fpe_t op)
{
    //TODO Check for XI
    fpe_set(rop->m_a, op);
    fpe_set(rop->m_b, op);
}

// Scalar multiple of an fp2e, store result in rop:
void fp2e_mul_fpe(fp2e_t rop, const fp2e_t op1, const fpe_t op2)
{
	fpe_mul(rop->m_a, op1->m_a, op2);
	fpe_mul(rop->m_b, op1->m_b, op2);
}

// Inverse multiple of an fp2e, store result in rop:
void fp2e_invert(fp2e_t rop, const fp2e_t op)
{
#ifdef BENCH
    numinvp2++;
    invp2cycles -= cpucycles();
#endif
	fpe_t tmp1, tmp2;  // Needed for intermediary results

	fpe_mul(tmp1, op->m_a, op->m_a);
	fpe_mul(tmp2, op->m_b, op->m_b);
#if(ALPHA == 2)
	fpe_double(tmp1, tmp1);
	fpe_sub(tmp2, tmp2, tmp1);
#elif(ALPHA == -2)
	fpe_double(tmp1, tmp1);
	fpe_add(tmp2, tmp2, tmp1);
#elif(ALPHA == -1)
	fpe_add(tmp2, tmp2, tmp1);
#else
#error "ALPHA must be -1, 2 or -2"
#endif
	fpe_invert(tmp2, tmp2);
	fpe_mul(rop->m_b, op->m_b, tmp2);
	fpe_neg(tmp2, tmp2);
	fpe_mul(rop->m_a, op->m_a, tmp2);
#ifdef BENCH
    invp2cycles += cpucycles();
#endif
}

// Print the fp2e:
void fp2e_print(FILE *outfile, const fp2e_t op)
{
	fprintf(outfile, "(");
	fpe_print(outfile, op->m_a);
	fprintf(outfile, " * a + \n");
	fpe_print(outfile, op->m_b);
	fprintf(outfile, ")");
}

typedef struct twistpoint_fp2_struct twistpoint_fp2_struct_t;

struct twistpoint_fp2_struct
{	
	fp2e_t m_x; // X-Coordinate (Jacobian Coordinate system)
	fp2e_t m_y; // Y-Coordinate (Jacobian Coordinate system)
	fp2e_t m_z; // Z-Coordinate (Jacobian Coordinate system)
	fp2e_t m_t; // T = Z^2, only used during pairing computation, set to zero if not set
};

typedef twistpoint_fp2_struct_t twistpoint_fp2_t[1];

curvepoint_fp_t curve_gen; // generator of E(\F_p)
twistpoint_fp2_t twist_gen; // generator of the subgroup of order n of E'(\F_{p^2}) 

void twistpoint_fp2_set(twistpoint_fp2_t rop, const twistpoint_fp2_t op)
{
	fp2e_set(rop->m_x, op->m_x);
	fp2e_set(rop->m_y, op->m_y);
	fp2e_set(rop->m_z, op->m_z);
  fp2e_setzero(rop->m_t);
}

void twistpoint_fp2_set_fp2e(twistpoint_fp2_t rop, const fp2e_t x, const fp2e_t y, const fp2e_t z)
{
	fp2e_set(rop->m_x, x);
	fp2e_set(rop->m_y, y);
	fp2e_set(rop->m_z, z);
  fp2e_setzero(rop->m_t);
}

void twistpoint_fp2_init_set_str(twistpoint_fp2_t rop, const char* xx, const char* xy, const char* yx, const char* yy){
	fp2e_t x,y,z;
	fp2e_set_str(x,xx,xy);
	fp2e_set_str(y,yx,yy);
	fp2e_setone(z);

	fp2e_set(rop->m_x, x);
	fp2e_set(rop->m_y, y);
	fp2e_set(rop->m_z, z);
    fp2e_setzero(rop->m_t);
}

void twistpoint_fp2_affineset_fp2e(twistpoint_fp2_t rop, const fp2e_t x, const fp2e_t y)
{
	fp2e_set(rop->m_x, x);
	fp2e_set(rop->m_y, y);
	fp2e_setone(rop->m_z);
  fp2e_setzero(rop->m_t);
}

void twistpoint_fp2_mixadd(twistpoint_fp2_t rop, const twistpoint_fp2_t op1, const twistpoint_fp2_t op2)
{
	fp2e_t tfp2e1, tfp2e2, tfp2e3, tfp2e4, tfp2e5, tfp2e6, tfp2e7, tfp2e8, tfp2e9; // Temporary variables needed for intermediary results
	fp2e_square(tfp2e1, op1->m_z);
	fp2e_mul(tfp2e2, op1->m_z, tfp2e1);
	fp2e_mul(tfp2e3, op2->m_x, tfp2e1);
	fp2e_mul(tfp2e4, op2->m_y, tfp2e2);
	fp2e_sub(tfp2e5, tfp2e3, op1->m_x);
	fp2e_sub(tfp2e6, tfp2e4, op1->m_y);
	fp2e_square(tfp2e7, tfp2e5);
	fp2e_mul(tfp2e8, tfp2e7, tfp2e5);
	fp2e_mul(tfp2e9, op1->m_x, tfp2e7);

	fp2e_double(tfp2e1, tfp2e9);
	fp2e_add(tfp2e1, tfp2e1, tfp2e8);
	fp2e_square(rop->m_x, tfp2e6);
	fp2e_sub(rop->m_x, rop->m_x, tfp2e1);
	fp2e_sub(tfp2e1, tfp2e9, rop->m_x);
	fp2e_mul(tfp2e2, tfp2e1, tfp2e6);
	fp2e_mul(tfp2e3, op1->m_y, tfp2e8);
	fp2e_sub(rop->m_y, tfp2e2, tfp2e3);
	fp2e_mul(rop->m_z, op1->m_z, tfp2e5);
}

void twistpoint_fp2_double(twistpoint_fp2_t rop, const twistpoint_fp2_t op)
{
	fp2e_t tfp2e1, tfp2e2, tfp2e3, tfp2e4; // Temporary variables needed for intermediary results
	fp2e_square(tfp2e1, op->m_y);
	fp2e_mul(tfp2e2, tfp2e1, op->m_x);
	fp2e_double(tfp2e2, tfp2e2);
	fp2e_double(tfp2e2, tfp2e2);
	fp2e_square(tfp2e3, tfp2e1);
	fp2e_double(tfp2e3, tfp2e3);
	fp2e_double(tfp2e3, tfp2e3);
	fp2e_double(tfp2e3, tfp2e3);
	fp2e_square(tfp2e4, op->m_x);
	fp2e_triple(tfp2e4, tfp2e4);
	fp2e_square(rop->m_x, tfp2e4);
	fp2e_double(tfp2e1, tfp2e2);
	fp2e_sub(rop->m_x, rop->m_x, tfp2e1);
	fp2e_sub(tfp2e1, tfp2e2, rop->m_x);
	fp2e_mul(rop->m_z, op->m_y, op->m_z);
	fp2e_double(rop->m_z, rop->m_z);
	fp2e_mul(rop->m_y, tfp2e4, tfp2e1);
	fp2e_sub(rop->m_y, rop->m_y, tfp2e3);
}

void twistpoint_fp2_mul(twistpoint_fp2_t rop, const twistpoint_fp2_t op, const mpz_t scalar)
{
	// TODO: Test...
	size_t i;
	twistpoint_fp2_t r;
	twistpoint_fp2_set(r, op);
	for(i = mpz_sizeinbase(scalar, 2) - 1; i > 0; i--)
	{
		twistpoint_fp2_double(r, r);
		if(mpz_tstbit(scalar, i - 1)) 
			twistpoint_fp2_mixadd(r, r, op);
	}
	twistpoint_fp2_set(rop, r);
}


void twistpoint_fp2_makeaffine(twistpoint_fp2_t op)
{
	fp2e_invert(op->m_z, op->m_z);
	fp2e_mul(op->m_y, op->m_y, op->m_z);
	fp2e_square(op->m_z, op->m_z);
	fp2e_mul(op->m_x, op->m_x, op->m_z);
	fp2e_mul(op->m_y, op->m_y, op->m_z);
	fp2e_setone(op->m_z);
}

void twistpoint_fp2_print(FILE *outfile, const twistpoint_fp2_t op)
{
	fprintf(outfile, "[");
	fp2e_print(outfile, op->m_x);
	fprintf(outfile, ", ");
	fp2e_print(outfile, op->m_y);
	fprintf(outfile, ", ");
	fp2e_print(outfile, op->m_z);
	fprintf(outfile, "]");
}

void points_init()
{
	fp2e_t twistgen_x;
	fp2e_t twistgen_y;
	fp2e_set_str(twistgen_x, BN_TWISTGEN_X);
	fp2e_set_str(twistgen_y, BN_TWISTGEN_Y);
	curvepoint_fp_set_str(curve_gen, BN_CURVEGEN);
	twistpoint_fp2_affineset_fp2e(twist_gen, twistgen_x, twistgen_y);
}

// Elements from F_{p^6}= F_{p^2}[Y] / (Y^3 - xi)F_{p^2}[Y] are represented as aY^2 + bY + c 
typedef struct fp6e_struct fp6e_struct_t;

struct fp6e_struct
{
	fp2e_t m_a;
	fp2e_t m_b;
	fp2e_t m_c;
};

typedef fp6e_struct_t fp6e_t[1];

fp2e_t xi; // constant coefficient in the irreducible polynomial Y^3 - xi, used to construct F_{p^6} as extension of F_{p^12}
fp2e_t ypminus1; // Y^{p-1} lies in F_{p^2}

fpe_t zeta; // Third root of unity in F_p fulfilling Z^{p^2} = -zeta * Z

// Multiples of powers of xi needed for cometa-pairing computation
fp2e_t xi2; // xi^2
fp2e_t _1o27xi3; // 1/27 xi^3
fp2e_t _1o3xi3; // 1/3 xi^3
fp2e_t _1o3xi; // 1/3 xi
fpe_t _1o3modp; // 1/3 \in \F_p

// Two constants needed for the cometa-pairing computation
fpe_t cometa_c0_const;
fpe_t cometa_c1_const;

void fp6_init()
{
	fp2e_set_str(xi, BN_XI);
	fp2e_set_str(ypminus1, BN_YPMINUS1);
	fpe_set_str(zeta, BN_ZETA);
	// fp2e_set_str(xi2, BN_XI2);
	// fp2e_set_str(_1o27xi3, BN_1O27XI3);
	// fp2e_set_str(_1o3xi3, BN_1O3XI3);
	// fp2e_set_str(_1o3xi, BN_1O3XI);
	// fpe_set_str(_1o3modp, BN_1O3MODP); // 1/3 \in \F_p
	
	// fpe_set_str(cometa_c0_const, BN_COMETA_C0_CONST);
	// fpe_set_str(cometa_c1_const, BN_COMETA_C1_CONST);
}

// Set fp6e_t rop to given value:
void fp6e_set(fp6e_t rop, const fp6e_t op)
{
	fp2e_set(rop->m_a, op->m_a);
	fp2e_set(rop->m_b, op->m_b);
	fp2e_set(rop->m_c, op->m_c);
}

// Initialize an fp6e, set to value given in three fp2es
void fp6e_set_fp2e(fp6e_t rop, const fp2e_t a, const fp2e_t b, const fp2e_t c)
{
	fp2e_set(rop->m_a, a);
	fp2e_set(rop->m_b, b);
	fp2e_set(rop->m_c, c);
}

// Initialize an fp6e, set to value given in six strings
void fp6e_set_str(fp6e_t rop, const char *a1, const char *a0, const char *b1, const char *b0, const char *c1, const char *c0)
{
	fp2e_set_str(rop->m_a, a1, a0);
	fp2e_set_str(rop->m_b, b1, b0);
	fp2e_set_str(rop->m_c, c1, c0);
}

// Set rop to one:
void fp6e_setone(fp6e_t rop)
{
	fp2e_setzero(rop->m_a);
	fp2e_setzero(rop->m_b);
	fp2e_setone(rop->m_c);
}

// Set rop to zero:
void fp6e_setzero(fp6e_t rop)
{
	fp2e_setzero(rop->m_a);
	fp2e_setzero(rop->m_b);
	fp2e_setzero(rop->m_c);
}

// Add two fp6e, store result in rop:
void fp6e_add(fp6e_t rop, const fp6e_t op1, const fp6e_t op2)
{
	fp2e_add(rop->m_a, op1->m_a, op2->m_a);
	fp2e_add(rop->m_b, op1->m_b, op2->m_b);
	fp2e_add(rop->m_c, op1->m_c, op2->m_c);
}

// Subtract op2 from op1, store result in rop:
void fp6e_sub(fp6e_t rop, const fp6e_t op1, const fp6e_t op2)
{
	fp2e_sub(rop->m_a, op1->m_a, op2->m_a);
	fp2e_sub(rop->m_b, op1->m_b, op2->m_b);
	fp2e_sub(rop->m_c, op1->m_c, op2->m_c);
}

// Subtract op2 from op1, store result in rop:
void fp6e_neg(fp6e_t rop, const fp6e_t op)
{
	fp2e_neg(rop->m_a, op->m_a);
	fp2e_neg(rop->m_b, op->m_b);
	fp2e_neg(rop->m_c, op->m_c);
}

// Multiply two fp6e, store result in rop:
void fp6e_mul(fp6e_t rop, const fp6e_t op1, const fp6e_t op2)
{
	fp2e_t tmp1, tmp2, tmp3, tmp4, tmp5, tmp6; // Needed for intermediary values

	// See "Multiplication and Squaring in Pairing-Friendly Fields", section 4, Karatsuba method
	fp2e_mul(tmp3, op1->m_a, op2->m_a);
	fp2e_mul(tmp2, op1->m_b, op2->m_b);
	fp2e_mul(tmp1, op1->m_c, op2->m_c);

	fp2e_add(tmp4, op1->m_a, op1->m_b);
	fp2e_add(tmp5, op2->m_a, op2->m_b);
	fp2e_mul(tmp6, tmp4, tmp5); 
	fp2e_sub(tmp6, tmp6, tmp2);
	fp2e_sub(tmp6, tmp6, tmp3);
	// fp2e_mulxi(tmp6, tmp6);
	fp2e_mul(tmp6, tmp6, xi);
	fp2e_add(tmp6, tmp6, tmp1);

	fp2e_add(tmp4, op1->m_b, op1->m_c);
	fp2e_add(tmp5, op2->m_b, op2->m_c);
	fp2e_mul(rop->m_b, tmp4, tmp5);
	fp2e_sub(rop->m_b, rop->m_b, tmp1);
	fp2e_sub(rop->m_b, rop->m_b, tmp2);
	// fp2e_mulxi(tmp4, tmp3);
	fp2e_mul(tmp4, tmp3, xi);
	fp2e_add(rop->m_b, rop->m_b, tmp4);

	fp2e_add(tmp4, op1->m_a, op1->m_c);
	fp2e_add(tmp5, op2->m_a, op2->m_c);

	fp2e_set(rop->m_c, tmp6);

	fp2e_mul(rop->m_a, tmp4, tmp5);
	fp2e_sub(rop->m_a, rop->m_a, tmp1);
	fp2e_add(rop->m_a, rop->m_a, tmp2);
	fp2e_sub(rop->m_a, rop->m_a, tmp3);
}

// Square an fp6e, store result in rop:
void fp6e_square(fp6e_t rop, const fp6e_t op)
{
	fp6e_mul(rop, op, op);
}
// Multiply with tau:
void fp6e_multau(fp6e_t rop, const fp6e_t op)
{
    fp2e_t tmp1;
    fp2e_set(tmp1, op->m_b);
    fp2e_set(rop->m_b, op->m_c);
    // fp2e_mulxi(rop->m_c, op->m_a);
	fp2e_mul(rop->m_c, op->m_a, xi);
    fp2e_set(rop->m_a, tmp1);
}

void fp6e_mul_fpe(fp6e_t rop, const fp6e_t op1, const fpe_t op2)
{
	fp2e_mul_fpe(rop->m_a, op1->m_a, op2);
	fp2e_mul_fpe(rop->m_b, op1->m_b, op2);
	fp2e_mul_fpe(rop->m_c, op1->m_c, op2);
}

void fp6e_mul_fp2e(fp6e_t rop, const fp6e_t op1, const fp2e_t op2)
{
	fp2e_mul(rop->m_a, op1->m_a, op2);
	fp2e_mul(rop->m_b, op1->m_b, op2);
	fp2e_mul(rop->m_c, op1->m_c, op2);
}

void fp6e_invert(fp6e_t rop, const fp6e_t op)
{
	fp2e_t tmp1, tmp2, tmp3, tmp4, tmp5;  // Needed to store intermediary results

	// See "Implementing cryptographic pairings"
	fp2e_square(tmp1, op->m_c);
	fp2e_mul(tmp5, op->m_a, op->m_b);
	// fp2e_mulxi(tmp5, tmp5);
	fp2e_mul(tmp5, tmp5, xi);
	fp2e_sub(tmp1, tmp1, tmp5); // A
	
	fp2e_square(tmp2, op->m_a);
	// fp2e_mulxi(tmp2, tmp2);
	fp2e_mul(tmp2, tmp2, xi);
	fp2e_mul(tmp5, op->m_b, op->m_c);
	fp2e_sub(tmp2, tmp2, tmp5); // B

	fp2e_square(tmp3, op->m_b);
	fp2e_mul(tmp5, op->m_a, op->m_c);
	fp2e_sub(tmp3, tmp3, tmp5); // C

	fp2e_mul(tmp4, tmp3, op->m_b);
	// fp2e_mulxi(tmp4, tmp4);
	fp2e_mul(tmp4, tmp4, xi);
	fp2e_mul(tmp5, tmp1, op->m_c);
	fp2e_add(tmp4, tmp4, tmp5);
	fp2e_mul(tmp5, tmp2, op->m_a);
	// fp2e_mulxi(tmp5, tmp5);
	fp2e_mul(tmp5, tmp5, xi);
	fp2e_add(tmp4, tmp4, tmp5); // F
	
	fp2e_invert(tmp4, tmp4);

	fp2e_mul(rop->m_a, tmp3, tmp4);
	fp2e_mul(rop->m_b, tmp2, tmp4);
	fp2e_mul(rop->m_c, tmp1, tmp4);
}

void fp6e_frobenius_p(fp6e_t rop, const fp6e_t op)
{
	fp2e_t tmp; // Needed to store intermediary results

	fp6e_set(rop, op);
	fpe_neg((rop->m_a)->m_a, (rop->m_a)->m_a);
	fpe_neg((rop->m_b)->m_a, (rop->m_b)->m_a);
	fpe_neg((rop->m_c)->m_a, (rop->m_c)->m_a);

	fp2e_mul(rop->m_b, rop->m_b, ypminus1);
	fp2e_square(tmp, ypminus1);
	fp2e_mul(rop->m_a, rop->m_a, tmp);
}

void fp6e_frobenius_p2(fp6e_t rop, const fp6e_t op)
{
	fpe_t tmp; // Needed for intermediary results

	fpe_square(tmp, zeta);
	fp2e_set(rop->m_c, op->m_c);
	fp2e_mul_fpe(rop->m_b, op->m_b, tmp);
	fp2e_mul_fpe(rop->m_a, op->m_a, zeta);
}

// Print the fp6e:
void fp6e_print(FILE * outfile, const fp6e_t op)
{
	// fprintf(outfile, "[");
	fp2e_print(outfile, op->m_a);
	fprintf(outfile, " * b^2 + \n");
	fp2e_print(outfile, op->m_b);
	fprintf(outfile, " * b + \n");
	fp2e_print(outfile, op->m_c);
	// fprintf(outfile, "]");
}

// Elements from F_{p^{12}}= F_{p^6}[Z] / (Z^2 - tau)F_{p^6}[Z] are represented as aZ + b
typedef struct fp12e_struct fp12e_struct_t;

struct fp12e_struct
{
	fp6e_t m_a;
	fp6e_t m_b;
};

typedef fp12e_struct_t fp12e_t[1];

fp6e_t tau; // F_{p^{12}} is constructed as F_{p^6}[Z]/(Z^2 - tau) F_{p^6}[Z]
fp2e_t zpminus1; // Z^{p-1}, lies in F_{p^2}
fp2e_t zpminus1inv; // Z^{p-1}, lies in F_{p^2}

void fp12_init()
{
	fp6e_set_str(tau, BN_TAU);
	fp2e_set_str(zpminus1, BN_ZPMINUS1);
	fp2e_set_str(zpminus1inv, BN_ZPMINUS1INV);
}

// Set fp12e_t rop to given value:
void fp12e_set(fp12e_t rop, const fp12e_t op)
{
	fp6e_set(rop->m_a, op->m_a);
	fp6e_set(rop->m_b, op->m_b);
}

// Initialize an fp12e, set to value given in two fp6es
void fp12e_set_fp6e(fp12e_t rop, const fp6e_t a, const fp6e_t b)
{
	fp6e_set(rop->m_a, a);
	fp6e_set(rop->m_b, b);
}

// Set rop to one:
void fp12e_setone(fp12e_t rop)
{
	fp6e_setzero(rop->m_a);
	fp6e_setone(rop->m_b);
}

// Set rop to zero:
void fp12e_setzero(fp12e_t rop)
{
	fp6e_setzero(rop->m_a);
	fp6e_setzero(rop->m_b);
}

// Add two fp12e, store result in rop:
void fp12e_add(fp12e_t rop, const fp12e_t op1, const fp12e_t op2)
{
	fp6e_add(rop->m_a, op1->m_a, op2->m_a);
	fp6e_add(rop->m_b, op1->m_b, op2->m_b);
}

// Subtract op2 from op1, store result in rop:
void fp12e_sub(fp12e_t rop, const fp12e_t op1, const fp12e_t op2)
{
	fp6e_sub(rop->m_a, op1->m_a, op2->m_a);
	fp6e_sub(rop->m_b, op1->m_b, op2->m_b);
}

// Multiply two fp12e, store result in rop:
void fp12e_mul(fp12e_t rop, const fp12e_t op1, const fp12e_t op2)
{
#ifdef BENCH
  nummultp12 ++;
  multp12cycles -= cpucycles();
#endif

	fp6e_t tmp1, tmp2, tmp3; // Needed to store intermediary results

	fp6e_mul(tmp1, op1->m_a, op2->m_a);
	fp6e_mul(tmp3, op1->m_b, op2->m_b);

	fp6e_add(tmp2, op2->m_a, op2->m_b);
	fp6e_add(rop->m_a, op1->m_a, op1->m_b);
	fp6e_set(rop->m_b, tmp3);

	fp6e_mul(rop->m_a, rop->m_a, tmp2);
	fp6e_sub(rop->m_a, rop->m_a, tmp1);
	fp6e_sub(rop->m_a, rop->m_a, rop->m_b);
	fp6e_multau(tmp1, tmp1);
	//fp6e_mul(tmp1, tmp1, tau);
	fp6e_add(rop->m_b, rop->m_b, tmp1);
#ifdef BENCH
  multp12cycles += cpucycles();
#endif
}

void fp12e_mul_fp6e(fp12e_t rop, const fp12e_t op1, const fp6e_t op2)
{
	fp6e_mul(rop->m_a, op1->m_a, op2);
	fp6e_mul(rop->m_b, op1->m_b, op2);
}

// Square an fp12e, store result in rop:
void fp12e_square(fp12e_t rop, const fp12e_t op)
{
#ifdef BENCH
  numsqp12 ++;
  sqp12cycles -= cpucycles();
#endif
	fp6e_t tmp1, tmp2, tmp3; // Needed to store intermediary results

	fp6e_mul(tmp1, op->m_a, op->m_b);

	fp6e_add(tmp2, op->m_a, op->m_b);
	fp6e_multau(tmp3, op->m_a);
	//fp6e_mul(tmp3, op->m_a, tau);
	fp6e_add(rop->m_b, tmp3, op->m_b);
	fp6e_mul(rop->m_b, rop->m_b, tmp2);

	fp6e_sub(rop->m_b, rop->m_b, tmp1);
	fp6e_multau(tmp2, tmp1);
	//fp6e_mul(tmp2, tmp1, tau);
	fp6e_sub(rop->m_b, rop->m_b, tmp2);

	fp6e_add(rop->m_a, tmp1, tmp1);
#ifdef BENCH
  sqp12cycles += cpucycles();
#endif
}

void fp12e_pow(fp12e_t rop, const fp12e_t op, const mpz_t exp)
{
	fp12e_t dummy;
	fp12e_set(dummy, op);
	fp12e_set(rop, op);
	int i;
	for(i = mpz_sizeinbase(exp, 2) - 1; i > 0; i--)
	{
		fp12e_square(rop, rop);
		if(mpz_tstbit(exp, i - 1)) 
			fp12e_mul(rop, rop, dummy);
	}
}


void fp12e_invert(fp12e_t rop, const fp12e_t op)
{
#ifdef BENCH
  numinvp12 ++;
  invp12cycles -= cpucycles();
#endif
	fp6e_t tmp1, tmp2; // Needed to store intermediary results

	fp6e_square(tmp1, op->m_a);
	fp6e_square(tmp2, op->m_b);
	fp6e_multau(tmp1, tmp1);
	//fp6e_mul(tmp1, tmp1, tau);
	fp6e_sub(tmp1, tmp2, tmp1);
	fp6e_invert(tmp1, tmp1);
	fp12e_set(rop, op);
	fp6e_neg(rop->m_a, rop->m_a);
	fp12e_mul_fp6e(rop, rop, tmp1);
#ifdef BENCH
  invp12cycles += cpucycles();
#endif
}

void fp12e_frobenius_p(fp12e_t rop, const fp12e_t op)
{
	fp6e_frobenius_p(rop->m_a, op->m_a);
	fp6e_frobenius_p(rop->m_b, op->m_b);
	fp6e_mul_fp2e(rop->m_a, rop->m_a, zpminus1);
}

void fp12e_frobenius_p2(fp12e_t rop, const fp12e_t op)
{
	fp6e_frobenius_p2(rop->m_a, op->m_a);
	fp6e_frobenius_p2(rop->m_b, op->m_b);
	fp6e_mul_fpe(rop->m_a, rop->m_a, zeta);
	fp6e_neg(rop->m_a, rop->m_a);
}

// Print the element to stdout:
void fp12e_print(FILE *outfile, const fp12e_t op)
{
	fp6e_print(outfile, op->m_a);
	fprintf(outfile, " * c + \n");
	fp6e_print(outfile, op->m_b);
}

void final_expo(fp12e_t rop)
{
	// First part: (p^6 - 1)
	fp12e_t dummy1, dummy2, dummy3, fp, fp2, fp3;
	fp12e_set(dummy1, rop);
	
	// This is exactly the p^6-Frobenius action:
	fp6e_neg(rop->m_a, rop->m_a);
	
	fp12e_invert(dummy2, dummy1);
	fp12e_mul(rop, rop, dummy2);

	// Second part: (p^2 + 1)
	fp12e_set(dummy1, rop);
	fp12e_frobenius_p2(rop, rop);
	fp12e_mul(rop, rop, dummy1);

	// Third part: Hard part (see Implementing cryptographic pairings over BN curves)
	fp12e_invert(dummy1, rop); // dummy1 = f^{-1}
	
	mpz_t exp;
	mpz_init_set(exp, x);
	mpz_mul_ui(exp, exp, 6);
	mpz_add_ui(exp, exp, 5);

	fp12e_pow(dummy2, dummy1, exp); // dummy2 = f^{-(6x+5)}
	
	fp12e_frobenius_p(dummy3, dummy2);
	fp12e_mul(dummy3, dummy3, dummy2); // dummy3 = f^{-(6x+5)p}*f^{-(6x+5)}
	fp12e_frobenius_p(fp, rop);
	fp12e_frobenius_p2(fp2, rop);
	fp12e_frobenius_p(fp3, fp2);
	
	fp12e_square(dummy1, rop);
	fp12e_square(dummy1, dummy1);

	fp12e_mul(rop, rop, fp); // rop = f*f^p

	mpz_set_ui(exp, 9);
	fp12e_pow(rop, rop, exp);
	fp12e_mul(rop, rop, dummy1);
	fp12e_mul(rop, rop, dummy2);
	fp12e_mul(rop, rop, dummy3);
	fp12e_mul(rop, rop, fp3);

	fp12e_square(dummy1, fp);
	fp12e_mul(dummy1, dummy1, fp2);
	fp12e_mul(dummy1, dummy1, dummy3);

	mpz_mul(exp, x, x);
	mpz_mul_ui(exp, exp, 6);
	mpz_add_ui(exp, exp, 1);
	fp12e_pow(dummy1, dummy1, exp);
	fp12e_mul(rop, rop, dummy1);

	mpz_clear(exp);
}

static void linefunction_add_ate(
        fp12e_t rop1, 
        twistpoint_fp2_t rop2, 
        const twistpoint_fp2_t op1, 
        const twistpoint_fp2_t op2, 
        const curvepoint_fp_t op3,
        const fp2e_t r2 /* r2 = y^2, see "Faster Computation of Tate Pairings" */
        )
{
    fp2e_t tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10; // Temporary variables needed for intermediary results
    fp6e_t tfp6e1, tfp6e2;


    fp2e_mul(tmp0, op2->m_x, op1->m_t); /* tmp0 = B = x2 * T1  = x2z1^2*/

    fp2e_add(tmp1, op2->m_y, op1->m_z);
    fp2e_square(tmp1, tmp1);
    fp2e_sub(tmp1, tmp1, r2);
    fp2e_sub(tmp1, tmp1, op1->m_t);
    fp2e_mul(tmp1, tmp1, op1->m_t); /* tmp1 = D = ((y2 + Z1)^2 - R2 - T1)T1  = 2y2z1^3 */

    fp2e_sub(tmp2, tmp0, op1->m_x); /* tmp2 = H = B - X1  = x2z1^2 - x1*/

    fp2e_square(tmp3, tmp2); /* tmp3 = I = H^2  = (x2z1^2 - x1)^2*/

    fp2e_double(tmp4, tmp3); 
    fp2e_double(tmp4, tmp4); /* tmp4 = E = 4I = 4(x2z1^2 - x1)^2*/

    fp2e_mul(tmp5, tmp2, tmp4); /* tmp5 = J = HE =  4(x2z1^2 - x1)(x2z1^2 - x1)^2*/

    fp2e_sub(tmp6, tmp1, op1->m_y); 
    fp2e_sub(tmp6, tmp6, op1->m_y); /* tmp6 = r = 2(D - 2Y1) = (2y2z1^3 - 2y1)*/
    
    fp2e_mul(tmp9, tmp6, op2->m_x); /* Needed later: tmp9 = x2(2y2z1^3 - 2y1)*/

    fp2e_mul(tmp7, op1->m_x, tmp4); /* tmp7 = V = X1*E = 4x1(x2z1^2 - x1)^2*/

    fp2e_square(rop2->m_x, tmp6);
    fp2e_sub(rop2->m_x, rop2->m_x, tmp5);
    fp2e_sub(rop2->m_x, rop2->m_x, tmp7);
    fp2e_sub(rop2->m_x, rop2->m_x, tmp7); /* X3 = r^2 - J - 2V = (2y2z1^3 - 2y1)^2 - 4(x2z1^2 - x1)(x2z1^2 - x1)^2 - 8x1(x2z1^2 - x1)^2*/

    fp2e_add(rop2->m_z, op1->m_z, tmp2);
    fp2e_square(rop2->m_z, rop2->m_z);
    fp2e_sub(rop2->m_z, rop2->m_z, op1->m_t);
    fp2e_sub(rop2->m_z, rop2->m_z, tmp3); /* Z3 = (z1 + H)^2 - T1 - I  = 2z1(x2z1^2 - x1) */
    
    fp2e_add(tmp10, op2->m_y, rop2->m_z); /* Needed later: tmp10 = y2 + z3*/

    fp2e_sub(tmp8, tmp7, rop2->m_x);
    fp2e_mul(tmp8, tmp8, tmp6);
    fp2e_mul(tmp0, op1->m_y, tmp5);
    fp2e_double(tmp0, tmp0);
    fp2e_sub(rop2->m_y, tmp8, tmp0); /* Y3 = r(V - X3) - 2Y1*J = (2y2z1^3 - 2y1)(4x1(x2z1^2 - x1)^2 - x3) - 8y1(x2z1^2 - x1)(x2z1^2 - x1)^2*/


    fp2e_square(rop2->m_t, rop2->m_z); /* T3 = Z3^2 */

    fp2e_square(tmp10, tmp10); /* tmp10 = (y2 + z3)^2 */
    fp2e_sub(tmp10, tmp10, r2);
    fp2e_sub(tmp10, tmp10, rop2->m_t); 
    fp2e_double(tmp9, tmp9);
    fp2e_sub(tmp9, tmp9, tmp10); /* tmp9 = 4x2(y2z1^3 - y1) - 2z3y2 */

    fp2e_mul_fpe(tmp10, rop2->m_z, op3->m_y); /* tmp10 = z3y_Q */
    fp2e_double(tmp10, tmp10);

    fp2e_neg(tmp6, tmp6);
    fp2e_mul_fpe(tmp1, tmp6, op3->m_x);
    fp2e_double(tmp1, tmp1);

    fp2e_setzero(tmp2);

    fp6e_set_fp2e(tfp6e1, tmp2, tmp9, tmp1);
    fp6e_set_fp2e(tfp6e2, tmp2, tmp2, tmp10);

    fp12e_set_fp6e(rop1, tfp6e1, tfp6e2);
}

static void linefunction_double_ate(fp12e_t rop1, twistpoint_fp2_t rop2, const twistpoint_fp2_t op1, const curvepoint_fp_t op3)
{
    fp2e_t tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp7, dummy; // Temporary variables needed for intermediary results
    fp6e_t tfp6e1, tfp6e2;

    fp2e_square(tmp0, op1->m_x); /* tmp0 = A = X1^2 = x1^2 */
    fp2e_square(tmp1, op1->m_y); /* tmp1 = B = Y1^2 = y1^2 */

    fp2e_square(tmp2, tmp1); /* tmp2 = C = B^2 = y1^4 */

    fp2e_add(tmp3, op1->m_x, tmp1);
    fp2e_square(tmp3, tmp3);
    fp2e_sub(tmp3, tmp3, tmp0);
    fp2e_sub(tmp3, tmp3, tmp2);
    fp2e_double(tmp3, tmp3); /* tmp3 = D = 2(X1 + B)^2 - A - C) = 4x1y1^2 */

    fp2e_triple(tmp4, tmp0); /* tmp4 = E = 3A = 3x1^2 */
    
    fp2e_add(tmp7, tmp4, op1->m_x); /* Needed later */

    fp2e_square(tmp5, tmp4); /* tmp5 = G = E^2 = 9x1^4 */

    fp2e_sub(rop2->m_x, tmp5, tmp3);
    fp2e_sub(rop2->m_x, rop2->m_x, tmp3); /* X3 = G - 2D = 9x1^4 - 8x1y1^2 */

    fp2e_add(rop2->m_z, op1->m_y, op1->m_z);
    fp2e_square(rop2->m_z, rop2->m_z);
    fp2e_sub(rop2->m_z, rop2->m_z, tmp1);
    fp2e_sub(rop2->m_z, rop2->m_z, op1->m_t); /* Z3 = (Y1 + Z1)^2 - B - T1 = 2y1z1; */

    fp2e_sub(rop2->m_y, tmp3, rop2->m_x);
    fp2e_mul(rop2->m_y, rop2->m_y, tmp4); 
    fp2e_double(dummy, tmp2);
    fp2e_double(dummy, dummy);
    fp2e_double(dummy, dummy);
    fp2e_sub(rop2->m_y, rop2->m_y, dummy); /* Y3 = E(D - X3) - 8C = 3x1^2(4x1y1^2 - X3) - 8y1^4 */

    fp2e_mul(tmp3, tmp4, op1->m_t);
    fp2e_double(tmp3, tmp3);
    fp2e_neg(tmp3, tmp3);
    fp2e_mul_fpe(tmp3, tmp3, op3->m_x); /* tmp3 = -6x1^2z1^2 * x_Q */

    fp2e_square(tmp7, tmp7);
    fp2e_sub(tmp7, tmp7, tmp0);
    fp2e_sub(tmp7, tmp7, tmp5);
    fp2e_double(dummy, tmp1);
    fp2e_double(dummy, dummy);
    fp2e_sub(tmp7, tmp7, dummy); /* tmp7 = 6x1^3 - 4y1^2 */

    fp2e_mul(tmp0, rop2->m_z, op1->m_t);
    fp2e_double(tmp0, tmp0);
    fp2e_mul_fpe(tmp0, tmp0, op3->m_y);
    
    fp2e_square(rop2->m_t, rop2->m_z); 

    fp2e_setzero(tmp1);

    fp6e_set_fp2e(tfp6e1, tmp1, tmp7, tmp3);
    fp6e_set_fp2e(tfp6e2, tmp1, tmp1, tmp0);

    fp12e_set_fp6e(rop1, tfp6e1, tfp6e2);
}

void optate(fp12e_t rop, const twistpoint_fp2_t op1, const curvepoint_fp_t op2)
{
    // op1 and op2 are assumed to be in affine coordinates!
    fpe_t zeta2;
    fp2e_t z2p, z3p;
    fpe_set_str(zeta2, BN_ZETA2);
    fp2e_set_str(z2p, BN_Z2P);
    fp2e_set_str(z3p, BN_Z3P);
    twistpoint_fp2_t q1, q2, q3;
    fp12e_setone(rop);


    fp12e_t dummy;
    fp2e_t tfp2e1, tfp2e2;

    mpz_t _6uplus2;
    mpz_init_set(_6uplus2, x);
    mpz_mul_ui(_6uplus2, _6uplus2, 6);
    mpz_add_ui(_6uplus2, _6uplus2, 2);

    twistpoint_fp2_t r, t;
    twistpoint_fp2_set(r, op1);
    fp2e_setone(r->m_t); /* As r has to be in affine coordinates this is ok */
    fp2e_setone(t->m_t); /* As t has to be in affine coordinates this is ok */

    fp2e_t r2;
    fp2e_square(r2, op1->m_y);

    unsigned long int i;

    unsigned long long int t1, t2, t3;
    // t1 = cpucycles();
    for(i = mpz_sizeinbase(_6uplus2, 2) - 1; i > 0; i--)
    {
        linefunction_double_ate(dummy, r, r, op2);
        fp12e_square(rop, rop);
        fp12e_mul(rop, rop, dummy);

        if (mpz_tstbit(_6uplus2, i - 1)) 
        {
            linefunction_add_ate(dummy, r, r, op1, op2, r2);
            fp12e_mul(rop, rop, dummy);
        }
    }
    // Needed because linfunction_add needs op2 in affine coordinates
    twistpoint_fp2_makeaffine(r);
    fp2e_setone(r->m_t); /* As r is in affine coordinates this is ok */

    fp2e_mul_fpe(tfp2e1, op1->m_x, zeta2);
    twistpoint_fp2_affineset_fp2e(q2, tfp2e1, op1->m_y); 

    fp2e_set(tfp2e1, op1->m_x);
    fpe_neg(tfp2e1->m_a, tfp2e1->m_a);
    fp2e_mul(tfp2e1, tfp2e1, z2p);
    fp2e_set(tfp2e2, op1->m_y);
    fpe_neg(tfp2e2->m_a, tfp2e2->m_a);
    fp2e_mul(tfp2e2, tfp2e2, z3p);
    twistpoint_fp2_affineset_fp2e(q1, tfp2e1, tfp2e2);

    fp2e_mul_fpe(tfp2e1, tfp2e1, zeta2);
    fp2e_neg(tfp2e2, tfp2e2);
    twistpoint_fp2_affineset_fp2e(q3, tfp2e1, tfp2e2);

    fp2e_setone(q3->m_t);
    fp2e_square(r2, q2->m_y);
    linefunction_add_ate(dummy, t, q3, q2, op2, r2);
    fp12e_mul(rop, rop, dummy);
    fp2e_square(r2, q1->m_y);
    linefunction_add_ate(dummy, t, t, q1, op2, r2);
    fp12e_mul(rop, rop, dummy);
    fp2e_square(r2, r->m_y);
    linefunction_add_ate(dummy, t, t, r, op2, r2);
    fp12e_mul(rop, rop, dummy);

    // t2 = cpucycles();
    final_expo(rop);
    // t3 = cpucycles();
    // printf("Miller: %llu\nFinal expo: %llu\n", t2-t1, t3-t2);
}

void ate(fp12e_t rop, const twistpoint_fp2_t op1, const curvepoint_fp_t op2)
{
    fp12e_setone(rop);

    fp12e_t dummy;

    mpz_t tminus1;
    mpz_init_set(tminus1, trace);
    mpz_sub_ui(tminus1, tminus1, 1);

    twistpoint_fp2_t r;
    twistpoint_fp2_set(r, op1);
    fp2e_setone(r->m_t); /* As r has to be in affine coordinates this is ok */
    
    fp2e_t r2;
    fp2e_square(r2, op1->m_y);

    unsigned long int i;

    // unsigned long long int t1, t2, t3;
    // t1 = cpucycles();
    for(i = mpz_sizeinbase(tminus1, 2) - 1; i > 0; i--)
    {
        linefunction_double_ate(dummy, r, r, op2);
        fp12e_square(rop, rop);
        fp12e_mul(rop, rop, dummy);

        if (mpz_tstbit(tminus1, i - 1)) 
        {
            linefunction_add_ate(dummy, r, r, op1, op2, r2);
            fp12e_mul(rop, rop, dummy);
        }
    }
    // t2 = cpucycles();
    final_expo(rop);
    // t3 = cpucycles();
    // printf("Miller: %llu\nFinal expo: %llu\n", t2-t1, t3-t2);
}

int fp12e_isone(fp12e_t op)
{	
	if(
	   !fpe_iszero(op->m_a->m_a->m_a) ||
	   !fpe_iszero(op->m_a->m_a->m_b) ||
	   !fpe_iszero(op->m_a->m_b->m_a) ||
	   !fpe_iszero(op->m_a->m_b->m_b) ||
	   !fpe_iszero(op->m_a->m_c->m_a) ||
	   !fpe_iszero(op->m_a->m_c->m_b) ||
	   !fpe_iszero(op->m_b->m_a->m_a) ||
	   !fpe_iszero(op->m_b->m_a->m_b) ||
	   !fpe_iszero(op->m_b->m_b->m_a) ||
	   !fpe_iszero(op->m_b->m_b->m_b) ||
	   !fpe_iszero(op->m_b->m_c->m_a) ||
	   !fpe_isone(op->m_b->m_c->m_b)
	   ){
		   return 0;
	   }
	   return 1;
}

int pairingProd2(curvepoint_fp_t a1, twistpoint_fp2_t b1, curvepoint_fp_t a2, twistpoint_fp2_t b2)
{
	fp12e_t res, e1, e2;
	optate(e1, b1, a1);
	optate(e2, b2, a2);
	fp12e_mul(res, e1, e2);

	return fp12e_isone(res);
}

int pairingProd3(curvepoint_fp_t a1, twistpoint_fp2_t b1, curvepoint_fp_t a2, twistpoint_fp2_t b2, curvepoint_fp_t a3, twistpoint_fp2_t b3)
{
	fp12e_t res, e1, e2, e3;
	optate(e1, b1, a1);
	optate(e2, b2, a2);
	optate(e3, b3, a3);
	fp12e_mul(res,e1,e2);
	fp12e_mul(res,res,e3);

	return fp12e_isone(res);
}

int main(int argc, char* argv[])
{
	// init
	fp_init();
	fp6_init();
	fp12_init();
	curve_init();
	points_init();

	// Input
	int num_elements = 18;
	int char_count = 79;
	int publicParamLeng = 2;
	mpz_t publicParam[publicParamLeng];
	curvepoint_fp_t proof_A, proof_Ap, proof_Bp, proof_C, proof_Cp, proof_H, proof_K;
	twistpoint_fp2_t proof_B;
	char **input = (char**) calloc(num_elements, sizeof(char*));
	for (int i=0;i< num_elements;i++)
	{
		input[i] = (char*) calloc(char_count, sizeof(char));
	}
	input[0] = "21139591401786686148758104085171563671377590446751309121396027748875156661644"; // "0x2ebc95b08158c9aec55a23806d77a17f69d540701ec56a318c8d76990bf1258c";
	input[1] = "16239878678216258978053794014066936183834024870847038813167529395815454742768";	// "0x23e77212cbec2e5b9f42d26c4990989e22a644ae4374a5b4d13a86b5a81e18f0";
	input[2] = "16965535099015174147383842414781227796161522935347727746754976911812233631503"; // "0x2582270f63c3a69f887c572411f5f47784e9de0918819e8d7ac4d9423ad9970f"; 
	input[3] = "15561464678002117391977951601264219454612210611207291050612306970299945975582"; // "0x22677a14f937a0dcda95069c1cfb622788489e30fa1c35ad528226e6862fff1e";
	input[4] = "5823845804124302971730045708813681134029160178246966568627155314774655374675"; // "0xce02e0ec5ebd85a0f4d6e1bc89092002d5ebe2d9be117972a199924d6a01553"; 
	input[5] = "1708374202911813861023252499661549051543719883493079056957775424728086391637"; // "0x3c6e7d102db86c82716a067cfba4b9850aa24ef0240528fd8ba12b03678ab55";
	input[6] = "7367543699215265691595216703788389867419984018028890936834266971345264404635"; // "0x1049e1c80bbe8b0c01b71c7cbcf5204455ba1cdec3ae74fd98006a270e79689b"; 
	input[7] = "11864092251170701746793868538144481285175488475786866257518327550479931988347"; // "0x1a3ad69780ef6e46a7e9f63019893c61162b37eef6f2a6ead9d73495f752d97b";
	input[8] = "14717677573969474141922286168032562895973221717491071827937179447310674073246"; // "0x2089e909cfe07a1c5dee80e56a5350a978a3dfc165c31409cae0d8e69384969e"; 
	input[9] = "20995607363248479762468353733132668173928517142046927728962546394297660669625"; // "0x2e6b17b7b119651b8cdab9159f5e829e64cbc3bff2ba42a6173c55914f24feb9";
	input[10] = "10278860691581916183915081202173847974947235523899090153863125732358434755061"; // "0x16b9a104fe24bc8908a5aef7650daa07f46e85e0142e2d2f8b2474295b0291f5"; 
	input[11] = "9034170296754417720730346712473312168011065464910915222054403525888321646886"; // "0x13f928c04f8ffe9520d425bdce6372423abc2e71b04a2886865a95c9cfedf126";
	input[12] = "12274159781618143788809767488697512140450633564234125255975902312714710303717"; // "0x1b22eda190a2a8593f78be31261a6290ddf656723bb5ab33292daa9818dc43e5"; 
	input[13] = "12123601003489349008509607354900740396982863236264158807606525155493969946366"; // "0x1acdb70a1b31d0c3e1a5751303edcf09fac5b45fed66ac8e1ee27bac550466fe";
	input[14] =	"11477826462049841305700340336476765269430941297673950876454374588555170341739"; // "0x19603835187801d3518e004570fee51ca17a42ecbe13e67e5c3799d6acb2e36b"; 
	input[15] = "6565678129488565228936863106582570930521016632442037468341748087033800287873"; // "0xe840ac90369827d2ba3a9a8b6f875de0cd31b2e552c1d0b5df337396071e681";
	input[16] = "18749801366771687119569815402027155847747319847751556037028760848497985844183"; // "0x297402f189b86081e344ce75cb0a089768156f1c72508d050c5650630265f3d7";
	input[17] = "9953815689031364393871846807411005549633476769935891842235399861038609338550"; // "0x1601a8f7c4010ac76e4d34416b364bf1667075fbee8c881d4eb38f0cd3e010b6";
	mpz_init_set_ui(publicParam[0], 113569);
	mpz_init_set_ui(publicParam[1], 1);

	curvepoint_fp_init_set_str(proof_A, input[0],input[1],"1");
	curvepoint_fp_init_set_str(proof_Ap, input[2],input[3],"1");
	twistpoint_fp2_init_set_str(proof_B,input[4],input[5],input[6],input[7]);
	curvepoint_fp_init_set_str(proof_Bp,input[8],input[9],"1");
	curvepoint_fp_init_set_str(proof_C,input[10],input[11],"1");
	curvepoint_fp_init_set_str(proof_Cp,input[12],input[13],"1");
	curvepoint_fp_init_set_str(proof_H, input[14],input[15],"1");
	curvepoint_fp_init_set_str(proof_K,input[16],input[17],"1");

	// verifying key
	curvepoint_fp_t vk_B, vk_gammaBeta1;
	twistpoint_fp2_t vk_A, vk_C, vk_gamma, vk_gammaBeta2, vk_Z;
	curvepoint_fp_t vk_IC[3];
	twistpoint_fp2_init_set_str(vk_A,"605277883741277067998757089620854213344342454131018764563731803821631833216","5780208356577781877029709601651611734344528471347057979924212196128087541173","3836234727767513625048875053821976424156523490917363464967806040709378810221","72616605417539107968777959647188823210945038483314474859581955894659533479");
	curvepoint_fp_init_set_str(vk_B,"5244050806582546246363707801080992827146630237438124204904608997543175263069","11043330697898197330582948926634797650426081926645217744080846758056305245009","1");
	twistpoint_fp2_init_set_str(vk_C,"19116878936370120509388482721749648506125663565436455480093376765454038283490","19853143919837539504459541822939260913892671116336907483088635058163949008109","5554626150521037890231274450708970331898920563100575922320903838566682594066","9744665906278717456709383416299220381585250396071523750881174529746981751802");
	twistpoint_fp2_init_set_str(vk_gamma,"4578504962507786833595020683331046674709651680203586493215769363853595471502","16442026267374291471227405674580212461254728302267323118746228740100681007950","12811615479952849803658431542199857970754267484025301652100963867304410263375","2125557773446690936767542125963571178846424772774481337624591517622295530193");
	curvepoint_fp_init_set_str(vk_gammaBeta1,"5704015748355616691766703021918762709403664425429126015526082821319544549091","21114552901101490257850821451640801607268410793266494132192179330860607294442","1");
	twistpoint_fp2_init_set_str(vk_gammaBeta2,"20356867301568381317194513063568441056883623538452547768726849395208639800758","18532479437056040881210070979066509018748404743746657283153535651012965996663","20136808714112366748646081419209125050695930963203151013425492891774449100640","1246889311156873502985742342422717621769073287813619502511915884369975675695");
	twistpoint_fp2_init_set_str(vk_Z,"2533779447106631414467935371715520351399707368570872000687283192664179011890","10829124150323656286722890006083703750633725979583862343624019688311264042421","5575488290461418253235413054264565713996725370582628939539397351274391964860","19273234910500420254672708417418248198395297039887155740386986936952927778252");
	curvepoint_fp_init_set_str(vk_IC[0],"7496470259527760717759918673269056479187764215826317640201523484569778378628","3129128929284658577384792948554307546849605941962627807241159293077937306055","1");
	curvepoint_fp_init_set_str(vk_IC[1],"18793153209555847362501657476662329963557481367201215580031331607748288985873","17996188574291334619261118489075220639857329901457072265162134727506284529000","1");
	curvepoint_fp_init_set_str(vk_IC[2],"7704112464825809742769231394063251018342023841535376926144202977714790216496","2918879334354371507852263189849384762743849561271180936966069170096867651260","1");

	// verify()
	curvepoint_fp_t vk_x, tmp0, tmp1, tmp2;
	curvepoint_fp_init_set_str(vk_x,"0","0","0");
	for(int i = 0; i< publicParamLeng; i++)
	{
		curvepoint_fp_mul(tmp0, vk_IC[i+1], publicParam[i]);
		curvepoint_fp_makeaffine(tmp0);
		if(fpe_iszero(vk_x->m_z)){
			curvepoint_fp_set(vk_x, tmp0);
		}
		else{
			curvepoint_fp_mixadd(vk_x, vk_x, tmp0);
			curvepoint_fp_makeaffine(vk_x);
		}
	}
	curvepoint_fp_mixadd(vk_x, vk_x, vk_IC[0]);
	curvepoint_fp_makeaffine(vk_x);

	curvepoint_fp_neg(tmp0, proof_Ap);
	curvepoint_fp_makeaffine(tmp0);
	if(!pairingProd2(proof_A, vk_A, tmp0, twist_gen)) return 1;
	
	curvepoint_fp_neg(tmp0, proof_Bp);
	curvepoint_fp_makeaffine(tmp0);
	if(!pairingProd2(vk_B, proof_B, tmp0, twist_gen)) return 2;
	
	curvepoint_fp_neg(tmp0, proof_Cp);
	curvepoint_fp_makeaffine(tmp0);
	if(!pairingProd2(proof_C,vk_C,tmp0, twist_gen)) return 3;
	
	curvepoint_fp_mixadd(tmp0, proof_A, proof_C);
	curvepoint_fp_makeaffine(tmp0);
	curvepoint_fp_mixadd(tmp0, vk_x, tmp0);
	curvepoint_fp_makeaffine(tmp0);
	curvepoint_fp_neg(tmp0, tmp0);
	curvepoint_fp_makeaffine(tmp0);
	curvepoint_fp_neg(tmp1, vk_gammaBeta1);
	curvepoint_fp_makeaffine(tmp1);
	if(!pairingProd3(proof_K, vk_gamma, tmp0, vk_gammaBeta2, tmp1, proof_B)) return 4;

	curvepoint_fp_mixadd(tmp0, vk_x, proof_A);
	curvepoint_fp_makeaffine(tmp0);
	curvepoint_fp_neg(tmp1, proof_H);
	curvepoint_fp_makeaffine(tmp1);
	curvepoint_fp_neg(tmp2, proof_C);
	curvepoint_fp_makeaffine(tmp2);
	if(!pairingProd3(tmp0, proof_B, tmp1, vk_Z, tmp2, twist_gen)) return 5;

	printf("transaction verified");

    return 0;
}