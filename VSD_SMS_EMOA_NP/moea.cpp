#include "algorithm.h"
#include <omp.h>

void InitializeBounds(int nvar, char * Instance)
{
	if( !strcmp("UF1", Instance) || !strcmp("UF2", Instance) || !strcmp("UF3", Instance) || !strcmp("UF4", Instance) || !strcmp("UF5", Instance) || !strcmp("UF6", Instance) || !strcmp("UF7", Instance) || !strcmp("UF8", Instance) || !strcmp("UF9", Instance) || !strcmp("UF10", Instance) || !strcmp("IMB1", Instance) || !strcmp("IMB2", Instance) || !strcmp("IMB3", Instance) || !strcmp("IMB4", Instance) || !strcmp("IMB5", Instance) || !strcmp("IMB6", Instance) || !strcmp("IMB7", Instance) || !strcmp("IMB8", Instance) || !strcmp("IMB9", Instance) || !strcmp("IMB10", Instance))
	{
		for(int i = 0 ;  i < nvar; i++)
		{
		   vlowBound[i]=0.0;
		   vuppBound[i]=1.0;//2.0*(i+1.0);
		}
	}
	
	if( !strcmp("WFG1", Instance) || !strcmp("WFG2", Instance) || !strcmp("WFG3", Instance) || !strcmp("WFG4", Instance) || !strcmp("WFG5", Instance) || !strcmp("WFG6", Instance) || !strcmp("WFG7", Instance) || !strcmp("WFG8", Instance) || !strcmp("WFG9", Instance))
	{
		for(int i = 0 ;  i < nvar; i++)
		{
		   vlowBound[i]=0.0;
		   vuppBound[i]=2.0*(i+1.0);
		}
	}
	if( !strcmp("DTLZ1", Instance) || !strcmp("DTLZ2", Instance) || !strcmp("DTLZ3", Instance) || !strcmp("DTLZ4", Instance) || !strcmp("DTLZ5", Instance) || !strcmp("DTLZ6", Instance) || !strcmp("DTLZ7", Instance) )
	{
		for(int i = 0 ;  i < nvar; i++)
		{
		   vlowBound[i]=0.0;
		   vuppBound[i]=1.0;
		}
	}
	if( !strcmp("RWP1", Instance))
        {
                for(int i = 0 ;  i < nvar; i++)
                {
                   vlowBound[i]=0.0;
                   vuppBound[i]=1.0;
                }
        }
        if( !strcmp("RWP2", Instance))
        {
                for(int i = 0 ;  i < nvar; i++)
                {
                   vlowBound[i]=1.0;
                   vuppBound[i]=3.0;
                }
        }

}
int main(int argc, char *argv[])
{

	int index = 1;
	int run = 1;
	strcpy(strpath, argv[index++]);
	strcpy(strTestInstance, argv[index++]);
	run= atoi(argv[index++]);
	nobj = atoi(argv[index++]);
	nPop= atoi(argv[index++]);
	nOffspring =  atoi(argv[index++]);
	max_nfes= atoll(argv[index++]);
	CR = atof(argv[index++]);
	realx=CR;
	F = atof(argv[index++]);
	nvar = atoi(argv[index++]);
        param_k =  (int)((4.0/24.0)*nvar);
        param_k = param_k - (int)(( nvar - param_k)%2);
        param_l = nvar-param_k;

	Di = sqrt(nvar)*atof(argv[index++]);
	Df = atof(argv[index++]);


	InitializeBounds(nvar, strTestInstance);

	CMOEAD MOEAD;
	MOEAD.exec_emo(run);
}
