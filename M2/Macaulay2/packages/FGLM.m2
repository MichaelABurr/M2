newPackage(
        "FGLM",
        Version => "1.0.0",
        Date => "January 29, 2019",
        Authors => {
	    { Name => "Dylan Peifer", Email => "djp282@cornell.edu", HomePage => "https://www.math.cornell.edu/~djp282" },
	    { Name => "Mahrud Sayrafi", Email => "mahrud@umn.edu", HomePage => "https://math.umn.edu/~mahrud" }
	    },
        Headline => "Compute Groebner bases via the FGLM algorithm",
        DebuggingMode => true
        )

-*
Copyright (C) 2019 Dylan Peifer and Mahrud Sayrafi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*-

export {"fglm"}

-------------------------------------------------------------------------------
--- top level functions
-------------------------------------------------------------------------------

---------------------------------------------------------------------------
-- See Section 2.4.4 of Thibaut Verron's thesis for details:
-- https://thibautverron.github.io/doc/2016-these.pdf
---------------------------------------------------------------------------
fglm = method()
fglm(Ideal,         Ring) := GroebnerBasis => (I1, R2) -> fglm(gb I1, R2)
fglm(GroebnerBasis, Ring) := GroebnerBasis => (G1, R2) -> (
    -- G1 = a Groebner basis for I1 in R1
    -- R2 = a polynomial ring
    -- returns Groebner basis G2 of the ideal generated by G1 in the ring R2

    R1 := ring G1;
    I1 := ideal gens G1; -- TODO: find a way to avoid recomputing GB
    if R1 == I1 then return forceGB(sub(gens G1, R2));
    if dim I1 > 0 then error "expected zero-dimensional ideal";
    if #gens R1 != #gens R2 then error "expected the same number of generators";

    M := multiplicationMatrices G1;
    m := numcols M#0;
    n := #gens R2;

    -- elements in the grobner basis
    G2 := {};
    -- elements in the staircase
    B2 := {1_R2};

    -- normal form translation table
    V := new MutableHashTable from {
	 B2#0 => transpose matrix { {1_(coefficientRing R1)} | toList ((m-1):0) } 
	 };

    -- list of elements between the staircase and grobner basis generators
    -- TODO: use a heap
    S := new MutableHashTable from for i to n - 1 list R2_i * 1_R2 => (i, 1_R2);

    while #S > 0 do (
	(elt, vals) := min pairs S;
	(i, mu) := vals;
	v := M#i * V#mu;

	VS := matrix { for s in B2 list V#s };
	if v % VS == 0 then (
--	if rank(VS | v) == numcols VS then ( -- not much faster
	    lambda := solve(VS, v);
--	    lambda := syz(VS | v, SyzygyRows => numcols VS); -- FIXME: why does this fail?
	    g := elt - matrix {B2} * lambda;
	    G2 = append(G2, g);
	    apply(pairs S, (k, v) -> if k % elt == 0 then remove(S, k));
	    ) else (
	    V#elt = v;
	    B2 = append(B2, elt);
	    for j to n - 1 do S#(R2_j * elt) = (j, elt);
	    apply(B2, b -> remove(S, b));
	    );
	);
    forceGB(matrix {G2})
    )

multiplicationMatrices = method()
multiplicationMatrices(GroebnerBasis) := List => (G) -> (
    -- G = a GroebnerBasis
    -- returns the matrices giving multiplication by variables in R/I

    R := ring G;
    I := ideal gens G;
    B := basis (R/I); -- TODO: find a way to avoid recomputing GB

    mons := first entries B;
    mats := {};

    for x in gens R do (
	F := first entries (x * B);
	M := transpose matrix for f in F list for m in mons list f_m;
	mats = append(mats, M);
    );

    mats
    )

cyclic = method(Options => {CoefficientRing => ZZ/32003, MonomialOrder => GRevLex})
cyclic(ZZ) := Ideal => opts -> (n) -> (
    R := (opts.CoefficientRing)[vars(0..n-1), MonomialOrder => opts.MonomialOrder];
    F := toList apply(1..n-1, d -> sum(0..n-1, i -> product(d, k -> R_((i+k)%n))))
         | {product gens R - 1};
    ideal F
    )

hcyclic = method(Options => {CoefficientRing => ZZ/32003, MonomialOrder => GRevLex})
hcyclic(ZZ) := Ideal => opts -> (n) -> (
    R := (opts.CoefficientRing)[vars(0..n), MonomialOrder => opts.MonomialOrder];
    F := toList apply(1..n-1, d -> sum(0..n-1, i -> product(d, k -> R_((i+k)%n)))) 
         | {product(n, i -> R_i) - R_n^n};
    ideal F
    )

katsura = method(Options => {CoefficientRing => ZZ/32003, MonomialOrder => GRevLex})
katsura(ZZ) := Ideal => opts -> (n) -> (
    n = n-1;
    R := (opts.CoefficientRing)[vars(0..n), MonomialOrder => opts.MonomialOrder];
    L := gens R;
    u := i -> (
	 if i < 0 then i = -i;
	 if i <= n then L_i else 0_R
	 );
    f1 := -1 + sum for i from -n to n list u i;
    F := toList prepend(f1, apply(0..n-1, i -> - u i + sum(-n..n, j -> (u j) * (u (i-j)))));
    ideal F
    )

-------------------------------------------------------------------------------
--- documentation
-------------------------------------------------------------------------------
beginDocumentation()

-------------------------------------------------------------------------------
--- tests
-------------------------------------------------------------------------------
test = (I1, MO2) -> (
    R1 := ring I1;
    R2 := (coefficientRing R1)(monoid ([gens R1], MonomialOrder => MO2));
    elapsedTime G2 := gb(sub(I1, R2));
    elapsedTime G2' := fglm(I1, R2);
    assert(gens G2 == gens G2')
    )

TEST ///
  debug needsPackage "FGLM"
  R1 = QQ[x,y,z]
  I1 = ideal(x*y + z - x*z, x^2 - z, 2*x^3 - x^2*y*z - 1)
  test(I1, Lex)
///

TEST ///
  debug needsPackage "FGLM"
  R1 = QQ[x,y,z]
  I1 = ideal(x^2 + 2*y^2 - y - 2*z, x^2 - 8*y^2 + 10*z - 1, x^2 - 7*y*z)
  test(I1, Lex)
///

TEST ///
  debug needsPackage "FGLM"
  R1 = QQ[x,y,z]
  I1 = ideal(x^2 + y^2 + z^2 - 2*x, x^3 - y*z - x, x - y + 2*z)
  test(I1, Lex)
///

TEST ///
  debug needsPackage "FGLM"
  R1 = QQ[x,y,z]
  I1 = ideal(x*y + z - x*z, x^2 - z, 2*x^3 - x^2*y*z - 1)
  test(I1, Lex)
///

TEST ///
  
///

end

restart
needsPackage "FGLM"
elapsedTime check FGLM -- ~2.7 seconds

restart
uninstallPackage "FGLM"
restart
installPackage "FGLM"
viewHelp "FGLM"

-- cyclic-6
-- gb: 0.310213
-- fglm: 4.19424
restart
needsPackage "FGLM"
load "examples.m2"
I = cyclic(6, MonomialOrder=>Lex)
G1 = elapsedTime gb I
I = cyclic(6)
R = newRing(ring I, MonomialOrder=>Lex)
G2 = elapsedTime fglm(I, R)


-- katsura-6
-- gb: 0.116438
-- fglm: 0.238552
restart
needsPackage "FGLM"
load "examples.m2"
I = katsura(6, MonomialOrder=>Lex)
G1 = elapsedTime gb I
I = katsura(6)
R = newRing(ring I, MonomialOrder=>Lex)
G2 = elapsedTime fglm(I, R)


-- cyclic-7
-- gb:
-- fglm:
restart
needsPackage "FGLM"
load "examples.m2"
I = cyclic(7, MonomialOrder=>Lex)
G1 = elapsedTime gb I
I = cyclic(7)
R = newRing(ring I, MonomialOrder=>Lex)
G2 = elapsedTime fglm(I, R)


-- katsura-7
-- gb: 6.82678
-- fglm: 1.252
restart
needsPackage "FGLM"
load "examples.m2"
I = katsura(7, MonomialOrder=>Lex)
G1 = elapsedTime gb I
I = katsura(7)
R = newRing(ring I, MonomialOrder=>Lex)
G2 = elapsedTime fglm(I, R)


-- katsura-8
-- gb: 
-- fglm: 
restart
needsPackage "FGLM"
load "examples.m2"
I = katsura(8, MonomialOrder=>Lex)
G1 = elapsedTime gb I
I = katsura(8)
R = newRing(ring I, MonomialOrder=>Lex)
G2 = elapsedTime fglm(I, R)


-- reimer-5
-- gb: 8.50152
-- fglm: 5.63204
restart
needsPackage "FGLM"
kk = ZZ/32003
R1 = kk[x,y,z,t,u, MonomialOrder=>Lex]
I1 = ideal(2*x^2 - 2*y^2 + 2*z^2 - 2*t^2 + 2*u^2 - 1,
           2*x^3 - 2*y^3 + 2*z^3 - 2*t^3 + 2*u^3 - 1,
           2*x^4 - 2*y^4 + 2*z^4 - 2*t^4 + 2*u^4 - 1,
           2*x^5 - 2*y^5 + 2*z^5 - 2*t^5 + 2*u^5 - 1,
           2*x^6 - 2*y^6 + 2*z^6 - 2*t^6 + 2*u^6 - 1)
G1 = elapsedTime gb I1
R2 = kk[x,y,z,t,u]
I2 = sub(I1, R2)
G2 = elapsedTime fglm(I2, R1)


-- virasoro
-- gb: 8.94464
-- fglm: 63.5117
restart
needsPackage "FGLM"
kk = ZZ/32003
R1 = kk[x1,x2,x3,x4,x5,x6,x7,x8, MonomialOrder=>Lex]
I1 = ideal(8*x1^2 + 8*x1*x2 + 8*x1*x3 + 2*x1*x4 + 2*x1*x5 + 2*x1*x6 + 2*x1*x7 - x1 - 8* x2*x3 - 2*x4*x7 - 2*x5*x6,
           8*x1*x2 - 8*x1*x3 + 8*x2^2 + 8*x2*x3 + 2*x2*x4 + 2*x2*x5 + 2*x2*x6 + 2*x2* x7 - x2 - 2*x4*x6 - 2*x5*x7,
	   -8*x1*x2 + 8*x1*x3 + 8*x2*x3 + 8*x3^2 + 2*x3*x4 + 2*x3*x5 + 2*x3*x6 + 2* x3*x7 - x3 - 2*x4*x5 - 2*x6*x7,
	   2*x1*x4 - 2*x1*x7 + 2*x2*x4 - 2*x2*x6 + 2*x3*x4 - 2*x3*x5 + 8*x4^2 + 8*x4* x5 + 2*x4*x6 + 2*x4*x7 + 6*x4*x8 - x4 - 6*x5*x8,
	   2*x1*x5 - 2*x1*x6 + 2*x2*x5 - 2*x2*x7 - 2*x3*x4 + 2*x3*x5 + 8*x4*x5 - 6*x4* x8 + 8*x5^2 + 2*x5*x6 + 2*x5*x7 + 6*x5*x8 - x5,
	   -2*x1*x5 + 2*x1*x6 - 2*x2*x4 + 2*x2*x6 + 2*x3*x6 - 2*x3*x7 + 2*x4*x6 + 2* x5*x6 + 8*x6^2 + 8*x6*x7 + 6*x6*x8 - x6 - 6*x7*x8,
	   -2*x1*x4 + 2*x1*x7 - 2*x2*x5 + 2*x2*x7 - 2*x3*x6 + 2*x3*x7 + 2*x4*x7 + 2* x5*x7 + 8*x6*x7 - 6*x6*x8 + 8*x7^2 + 6*x7*x8 - x7,
	   -6*x4*x5 + 6*x4*x8 + 6*x5*x8 - 6*x6*x7 + 6*x6*x8 + 6*x7*x8 + 8*x8^2 - x8)
G1 = elapsedTime gb I1
R2 = kk[x1,x2,x3,x4,x5,x6,x7,x8]
I2 = sub(I1, R2)
G2 = elapsedTime fglm(I2, R1)


-- chemkin
-- gb: 
-- fglm: 
restart
needsPackage "FGLM"
kk = ZZ/32003
R1 = kk[w,x3,x4,y2,y3,y4,y5,z2,z3,z4,z5, MonomialOrder=>Lex]
I1 = ideal(-4*w*y2 + 9*y2^2 + z2,
           x3^2 + y3^2 + z3^2 - 1,
           x4^2 + y4^2 + z4^2 - 1,
           9*y5^2 + 9*z5^2 - 8,
           -6*w*x3*y2 + 3*x3 + 3*y2*y3 + 3*z2*z3 - 1,
           3*x3*x4 + 3*y3*y4 + 3*z3*z4 - 1,
           x4 + 3*y4*y5 + 3*z4*z5 - 1,
           -6*w + 3*x3 + 3*x4 + 8,
           9*y2 + 9*y3 + 9*y4 + 9*y5 + 8,
           z2 + z3 + z4 + z5,
           w^2 - 2)
G1 = elapsedTime gb I1
R2 = kk[x1,x2,x3,x4,x5,x6,x7,x8]
I2 = sub(I1, R2)
G2 = elapsedTime fglm(I2, R1)
