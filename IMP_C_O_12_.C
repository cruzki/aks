//Librería con la implementación del AKS (test de primalidad)

#include <NTL/ZZXFactoring.h> //libreria para factorizar en Z_p[x]
#include <NTL/ZZ_pEX.h> //libreria para computar en Z_p[x]/x^r-1

NTL_CLIENT

long orden_G(ZZ n, long r)
//funcion que calcula el orden de n modulo r, o_r(n)
// !!!!!!! NO CONTROLA LA EXISTENCIA !!!!!!!!
//if (GCD(n,to_ZZ(r))!=1) ord=0; //compruebo si tiene orden
{
//declaración de variables

	ZZ_p::init(to_ZZ(r)); //inicializo el modulo ZZ_p
	ZZ_p aux = to_ZZ_p(n); //auxiliar para ir calculando potencias inicializado a n 
    long ord; //auxiliar para calcular el orden
 
   ord=1; //inicializo
   while (aux != 1)
   {
      mul(aux,aux,to_ZZ_p(n)); //aux = aux*n  mod r = n^ord
      ord++;
   }
   return (ord);
}

int AKS(ZZ n, int opcion)
//función que implementa un test determínistico de primalidad AKS.
//si opcion es igual a 0 se realiza el test estandar, para cualquier
//otro valor usamos la conjetura propuesta por AKS que hace el 
//algoritmo O(log^3(n)).
//Devuelve 0 en caso de ser compuesto, 1 si es primo.
//Actualmente sólo podemos trabajar hasta entero de tamaño long
//por problemas con el exponente :(
{
//declaracion de variables

	unsigned long i,j; //contadores, pueden llegar a ser muy grande
	unsigned long r; //módulo sobre el que dividimos
	unsigned long ord; //auxiliar para guardar los ordenes 
    ZZ aux; //auxiliar para cálculos en la parte de la conjetura
	ZZX f; //polinomio a factorizar en paso 1.
	ZZ_pX modulo; //modulo sobre el que trabajaremos en el paso 5
	ZZ_pX pot1,pot2; //auxiliar para calcular las potencias del paso 5
	ZZ_pXModulus Modulo; //variable en la que guardaremos precalculos
	vec_pair_ZZX_long factores; //factorización, me la devuelve asi

   	//Paso 1

   i=2;
   while (i < trunc(log(n))+1)
   {
      f=ZZX(i,1) - n; //polinomio para hayar la factorización, x^i-n
      ZZ c; //variable que contendra el mcd de los coeficiente... en este caso basura     
      factor(c,factores,f,0,0); //factorizo
      for (j = 1; j < factores.length()+1; j++) 
      {
         if ((deg(factores(j).a) == 1) & (factores(j).b != i)) return 0;
      } //compruebo grado de los factores
      i++;
   } //si sale de aki no existe a,i tales que a^i=n.

	//Paso 2
   if (opcion==0)
   {
      r=2;
      if (GCD(n,to_ZZ(r))!=1) ord=0; //compruebo si tiene orden
      else ord=orden_G(n,r);  //calulo
      while (ord < trunc(4*log(n)*log(n))+1)
      {
         r++;
         if (GCD(n,to_ZZ(r))!=1) ord=0;
         else ord=orden_G(n,r);
      } //obtengo r tal que o_r(n)>4log^2(n)
   } //algoritmo clásico.
   else
   {
      r=2;
      sqr(aux,n); //aux=n^2
      sub(aux,n,1); //aux=n^2-1
      while (divide(aux,r)) r++; // r=n^2-1?
   }//usando conjetura para O(log^3(n))
    
	//Paso 3
   for (i=2; i < r; i++) if ((GCD(to_ZZ(i),n)!=1) & (i < n))  return 0;
     
  	//Paso 4
   
   if (n < r) return 1;

  	//Paso 5
   
   ZZ_p::init(n); //inicializo el modulo Z_p
   i=1;
   modulo=ZZ_pX(r, to_ZZ_p(1)); //x^r
   sub(modulo,modulo,1); //x^r-1
   build(Modulo,modulo); //construyo precalculos para x^r-1
   
   while (i < 2*SqrRoot(r)*log(n))
   {
      PowerXPlusAMod(pot1,to_ZZ_p(i),n,Modulo); // (x+i)^n mod x^r-1,n
      PowerXMod(pot2,n,Modulo); //x^n mod x^r-1,n
      if (!(IsZero((pot1 - pot2 - i) % Modulo))) return(0);
      i++; //compruebo (x+i)^n-x^n-i=0 mod x^r-1,n 
   }   
	//Paso 6
	
   return(1);
}  
