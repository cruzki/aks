#include <cstdlib>
#include <iostream>
#include <NTL/ZZ.h> //librería con números ZZ
#include "imp_C_O(12).h" //librería con el AKS
#include <fstream> //escritura de archivos
NTL_CLIENT

int main()
//Programa que realiza una prueba sobre la mejora del
//la probabilidad de error del algoritmo de Miller-Rabin
//al aumentar el numero de iteraciones.
{
    ZZ a; //inicio
    ZZ i; //contador
    ZZ b; //final
    unsigned long j; //contador de primos
    unsigned long k; //contador
    short aux; //auxiliar
    
    ifstream entrada("rango.txt", ios::in); //abro fichero con rango 
    entrada >> a >> b; //leo limites
    entrada.close();   //cierro fichero
    
    ofstream salida;   // objeto de la clase ofstream
    salida.open("salida.txt",ios::out); //abro archivo

for (k=1; k <= b; k++)
{
   j=0;
   for (i=1; i <=100000; i++)  //criba probabilística
   {  
      aux=ProbPrime(to_ZZ(a), k);
      j+=aux;
   } 
   
   salida << k << " " << j << endl;
}    
   salida.close(); //cierro el fichero
   
} // main
