// Cabezeras de la libreria con el algoritmo de primalidad AKS

#include <NTL/ZZ_pXFactoring.h> //libreria para factorizar en Z_p[x]
#include <NTL/ZZ_pEX.h> //libreria para computar en Z_p[x]/x^r-1

NTL_CLIENT

ZZ orden(ZZ n, ZZ r); //calcula el orden de n modulo r, o_r(n)
int AKS(ZZ n, int opcion); //test de primalidad AKS
