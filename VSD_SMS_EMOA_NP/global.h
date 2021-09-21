#ifndef __GLOBAL_H_
#define __GLOBAL_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
#include <cassert>
#include <algorithm>


using namespace std;

#include "random.h"


//------------- Parameters in test instance ------------------
enum DominanceRelation
{
        INCOMPARABLE = 0,       // LHS and RHS are incomparable
        LHS_DOMINATES_RHS = 1,  // LHS strictly dominates RHS
        RHS_DOMINATES_LHS = 2,  // RHS strictly dominates LHS
        EQUIVALENT = 3,         // LHS and RHS are equally good
};


int     nvar,  nobj;                    //  the number of variables and objectives
int param_k, param_l;

double  lowBound = 0,   uppBound = 1;   //  lower and upper bounds of variables
double  vlowBound[2000] ,   vuppBound[2000];   //  lower and upper bounds of variables

char    strTestInstance[256];
char    strpath[800];


//------------- Parameters in random number ------------------
int     seed    = 177;
long    rnd_uni_init;        


//------------- Parameters in MOEA/D -------------------------

double          scale[100];  


int		etax    = 2, 	etam    = 50;   // distribution indexes of crossover and mutation

double  realx,  realm,  realb = 0.9;    // crossover, mutation, selection probabilities

double Di, Df, CR, F; // distance available of the hypersphere...
int nPop, nOffspring;
long long max_nfes;
inline DominanceRelation dominance(vector<double> const &lhs, vector<double> const &rhs){
   std::size_t l = 0, r = 0;
   int nobj=lhs.size();
        for (std::size_t i=0; i< nobj; i++)
        {
                if (lhs[i] < rhs[i]) l++;
                else if (lhs[i] > rhs[i]) r++;
        }

        if (l > 0)
        {
                if (r > 0) return INCOMPARABLE;
                else return LHS_DOMINATES_RHS;
        }
        else
        {
                if (r > 0) return RHS_DOMINATES_LHS;
                else return EQUIVALENT;
        }

}

#endif
