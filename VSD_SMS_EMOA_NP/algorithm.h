#ifndef __EVOLUTION_H_
#define __EVOLUTION_H_

#include <cfloat>
#include <set>
#include <queue>
#include <map>
#include <unordered_set>
#include <iomanip>
#include "global.h"
#include "recomb.h"
#include "individual.h"
#include "HypervolumeIndicator.h"
using namespace shark;
class CMOEAD
{

   public:
	CMOEAD();
	virtual ~CMOEAD();

	void init_population();                
	void replacement_phase();
	void evol_population();                                    
	// execute MOEAD
	void exec_emo(int run);
	void save_front(char savefilename[4024]); 
	void save_pos(char savefilename[4024]);
	void update_parameterD();
	void reference_point();
	double distance_var( vector<double> &a, vector<double> &b);
	void full_dominance_information();
	void classic_hv_selection_diversity();
	void update_archive();
   private:
	vector <CIndividual> pool;
	vector<vector<double> >archive_yobj, archive_xvar;
        vector<int> Np,Rp;//rank
	vector<unordered_set<int> > Sp;//dominated indexes and inverse
        vector<vector<int> > fronts ;
  	HypervolumeIndicator m_indicator, m_indicator2;

	// algorithm parameters
	long long nfes;          //  the number of function evluations
	double	D;	//Current minimum distance

};
CMOEAD::CMOEAD()
{

}
CMOEAD::~CMOEAD()
{

}
void CMOEAD::update_parameterD()
{
      double TElapsed = nfes, TEnd = max_nfes;
      D = Di - Di * (TElapsed / (TEnd*Df));
      //D=max(D, 0.0);
}
double CMOEAD::distance_var( vector<double> &a, vector<double> &b)
{
   double dist = 0 ;
   for(int i = 0; i < a.size(); i++)
   {
      double factor = (a[i]-b[i])/(vuppBound[i]-vlowBound[i]);
      dist += factor*factor;
   }
   return sqrt(dist);
}
void CMOEAD::init_population()
{
     vector<double> reference(nobj, DBL_MAX);
     m_indicator.setReference(reference);
     m_indicator2.setReference(reference);
    for(int i=0; i< nPop; i++)
    {
	       CIndividual ind;
		// Randomize and evaluate solution
		ind.rnd_init();
		ind.obj_eval();
		// Initialize the reference point
		pool.push_back(ind); 
		archive_yobj.push_back(ind.y_obj);
		archive_xvar.push_back(ind.x_var);
		//   nfes++;
     }
    reference_point();
}
void CMOEAD::evol_population()
{
//	static bool renovate =false;
//	if(D<0 && !renovate){
//		for(int i = 0; i < nPop; i++){
//		   pool[i].y_obj=archive_yobj[i];
//		   pool[i].x_var=archive_xvar[i];
//		}
//		renovate=true;
//	}
   full_dominance_information();
       int idx1=rand()%nPop, idx2=rand()%nPop, idx3=rand()%nPop, idx4=rand()%nPop;
      if(Rp[idx2] < Rp[idx1]) idx1=idx2;
      else if(Rp[idx1] == Rp[idx2]) idx1=(rand()%2)?idx1:idx2;

      if(Rp[idx4] < Rp[idx3]) idx3=idx4;
      else if(Rp[idx3] == Rp[idx4]) idx3=(rand()%2)?idx3:idx4;
      // produce a child solution
      CIndividual child1 = pool[idx1], child2 = pool[idx3];

      bool crosed = real_sbx_xoverA(pool[idx1], pool[idx3], child1, child2);
      // apply polynomial mutation
      double pm=rnd_uni(&rnd_uni_init);
    if (pm <= 0.5){
        bool mut1 = realmutation(child1, 1.0/(double)nvar);
        child1.obj_eval();
	if(crosed || mut1)nfes++;
	pool.push_back(child1);
	archive_yobj.push_back(child1.y_obj);
	archive_xvar.push_back(child1.x_var);
    }
    else{
        bool mut2 = realmutation(child2, 1.0/(double)nvar);
        child2.obj_eval();
	if(crosed || mut2)nfes++;
	pool.push_back(child2);
	archive_yobj.push_back(child2.y_obj);
	archive_xvar.push_back(child2.x_var);
   }
   reference_point();
   full_dominance_information();
   classic_hv_selection_diversity();
   update_archive();
}
void CMOEAD::exec_emo(int run)
{
        char filename1[5024];
        char filename2[5024];
	seed = run;
	srand(run);
	seed = (seed + 23)%1377;
	rnd_uni_init = -(long)seed;

	//initialization
	nfes      = 0;
	init_population();

	sprintf(filename1,"%s/POS/POS_R2_EMOA_%s_RUN%d_seed_%d_nobj_%d_nvar_%d_DI_%lf_DF_%lf_Px_%lf",strpath, strTestInstance,run, seed, nobj, nvar, Di/sqrt(nvar), Df, CR);
	sprintf(filename2,"%s/POF/POF_R2_EMOA_%s_RUN%d_seed_%d_nobj_%d_nvar_%d_DI_%lf_DF_%lf_Px_%lf",strpath, strTestInstance,run, seed, nobj, nvar, Di/sqrt(nvar), Df, CR);
        long long current = nfes;
	long long accumulator = 0, bef = nfes;
	save_pos(filename1);
        save_front(filename2);
	while(nfes<max_nfes)
	{
		update_parameterD();
		evol_population();
		accumulator += nfes - bef ;
                if(accumulator > 0.1*(max_nfes)  )
		{
	           accumulator -= 0.1*(max_nfes);
		   save_pos(filename1);
		   save_front(filename2);
		}
		bef=nfes;
	        //nfes += nOffspring;
	}
	save_pos(filename1);
	save_front(filename2);
}
void CMOEAD::save_front(char saveFilename[4024])
{
    std::fstream fout;
    fout.open(saveFilename,fstream::app|fstream::out );
    for(int n=0; n < nPop; n++)
    {
       for(int k=0;k<nobj;k++)
          fout<<archive_yobj[n][k]<<"  ";
       for(int k=0;k<nobj;k++)
          fout<<m_indicator.m_reference[k]<<"  ";
       for(int k=0;k<nobj;k++)
          fout<<m_indicator2.m_reference[k]<<"  ";

       for(int k=0;k<nobj;k++)
          fout<<pool[n].y_obj[k]<<"  ";
       fout<<"\n";
    }
    fout.close();
}

void CMOEAD::save_pos(char saveFilename[4024])
{
   std::fstream fout; //fout.open(saveFilename,std::ios::out);
   fout.open(saveFilename, fstream::app|fstream::out);
   for(int n=0; n<nPop; n++)
   {

      for(int k=0;k<nvar;k++)
         fout<<archive_xvar[n][k] << "  ";

      for(int k=0;k<nvar;k++)
         fout<<pool[n].x_var[k] << "  ";

   //   for(int k=0;k<nvar;k++)
   //	 fout<<pool[parent_idx[n]].x_var[k]/(vuppBound[k]-vlowBound[k]) << "  "; //fout<<population[n].indiv.x_var[k]<< fixed << setprecision(30) << "  ";
   	fout<<"\n";
   }
   fout.close();
}
void CMOEAD::reference_point(){
  vector<double> reference(nobj, 0);
  for(int m =0; m < nobj; m++){
    for(auto &ind:pool){
	   reference[m]=max(reference[m], ind.y_obj[m]);
    }
    reference[m] +=1.0;
  }
  m_indicator.setReference(reference);
}
void CMOEAD::full_dominance_information()
{
   int n=pool.size();
   Sp.assign(n, unordered_set<int>());
   Np.assign(n, 0);
   Rp.assign(n, 0);
   fronts.assign(1, vector<int>());
   int rank = 0;
   for(int pidx1=0; pidx1< n; pidx1++)
   {
      for(int pidx2=0; pidx2< n; pidx2++)
      {
	if(pidx1 == pidx2) continue;
	DominanceRelation ref = dominance(pool[pidx1].y_obj, pool[pidx2].y_obj);
        if(ref==LHS_DOMINATES_RHS) Sp[pidx1].insert(pidx2);
 	else if(ref==RHS_DOMINATES_LHS) Np[pidx1]++;

 //       if( pool[pidx1] < pool[pidx2]) Sp[pidx1].insert(pidx2);
 //	else if( pool[pidx2] < pool[pidx1]) Np[pidx1]++;
      }
      if( Np[pidx1] == 0)
      {
         fronts[rank].push_back(pidx1);
	 Rp[pidx1]=rank;
      }
   }
   while(true)
   {
      vector<int> next_front;
      for(auto idx:fronts[rank])
      {
	for(auto idx_dominated:Sp[idx])
        {
	  Np[idx_dominated]--;
          if(Np[idx_dominated]  == 0) 
          {
	     next_front.push_back(idx_dominated);
	     Rp[idx_dominated] = rank+1;
          }
        }
      }
      if(next_front.empty()) break;
      fronts.push_back(next_front);
      rank++;
   }
}
void CMOEAD::classic_hv_selection_diversity()
{
  if(D>0){
     //find the nearest pair...
     pair<double, pair<int, int> >pp(DBL_MAX, make_pair(-1, -1));
     for(int i = 0; i < pool.size(); i++){
       for(int j = i+1; j < pool.size(); j++){
   	double dist = distance_var(pool[i].x_var, pool[j].x_var);
   	pp=min(pp, make_pair(dist, make_pair(i,j)));
       }
     }
     if(pp.first<=D){
       int idx_to_remove=-1;
       if( Rp[pp.second.first] > Rp[pp.second.second] )idx_to_remove=pp.second.first;
       else if( Rp[pp.second.first] < Rp[pp.second.second] )idx_to_remove=pp.second.second;
       else{
   	 int rank = Rp[pp.second.first];
	 sort(fronts[rank].begin(), fronts[rank].end());
            vector<vector<double> > lastfront;
            for(auto idx:fronts[rank])lastfront.push_back(pool[idx].y_obj);
       	 vector<pair<double, size_t> > to_remove=m_indicator.leastContributors(lastfront, fronts[rank].size());
   	 for(auto i:to_remove){
   		 if(fronts[rank][i.second]==pp.second.first){idx_to_remove=pp.second.first; break;}
   		 if(fronts[rank][i.second]==pp.second.second){ idx_to_remove=pp.second.second; break;}
   	 }
       }
       pool[idx_to_remove] = pool.back();
       pool.pop_back();
       return;
     }
 }

  vector<vector<double> > lastfront;
   sort(fronts.back().begin(), fronts.back().end());
  for(auto idx:fronts.back()) lastfront.push_back(pool[idx].y_obj);
    vector<pair<double, size_t> > to_remove=m_indicator.leastContributors(lastfront, 1);
    int idx_to_remove=fronts.back()[to_remove[0].second];
    pool[idx_to_remove] = pool.back();
    pool.pop_back();
}
void CMOEAD::update_archive(){
	if(archive_yobj.size()<nPop)return;
  vector<double> reference(nobj, 0);
  for(int m =0; m < nobj; m++){
    for(auto &point:archive_yobj){
	   reference[m]=max(reference[m], point[m]);
    }
    reference[m] +=1.0;
  }
  m_indicator2.setReference(reference);

    int idx=m_indicator2.leastContributors(archive_yobj, 1)[0].second;
    archive_yobj[idx]=archive_yobj.back();
    archive_yobj.pop_back();
    archive_xvar[idx]=archive_xvar.back();
    archive_xvar.pop_back();
}
#endif
