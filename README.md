# chem-eng
Joint project with Chem Eng
Initial C++ Tabu Code Files:
Possibly imcomplete, I grabbed these two files friday afternoon prior to having a full understanding of the program
If I am able to determine that more files are needed to run the optimization algorithm, I willl track them down and add them
-Peter

tabusearch.cpp
#include "pd2_tabuSearch.hpp"
#include "cmath"

pd2_tabuSearch::pd2_tabuSearch(){
  neighbor_list = NULL;
  tabu_list = NULL;
  cur_sol = NULL;
  init_sol = NULL;
  best_sol = NULL;
  nni = 0;
  nni_limit = 0;
  calculateObj = NULL;
  cmpTabu = NULL;
  generateNeighbor = NULL;
}

double pd2_tabuSearch::initialObjective(){
  calculateObj( init_sol );
  return ( init_sol->obj );
}

void pd2_tabuSearch::search(){
  solution_t *tmpSolPtr; // for swapping solutions
  std::list<solution_t*>::iterator nl_it;
  std::list<solution_t*>::iterator tl_it;
  bool isTabu;

  calculateObj(init_sol);
  //std::cout << "Initial objective function: " << init_sol->obj << std::endl;
  assignSol(init_sol, cur_sol);
  assignSol(init_sol, best_sol);
  //just assign the initial solution to the whole tabu list
  //that will slow down the first sevral iterations a little
  //checking the same thing needlessly, but it saves me from
  //having to worry if the list is full or not.
  for(tl_it = tabu_list->begin(); tl_it != tabu_list->end(); ++tl_it){
    assignSol(init_sol, *tl_it);
  }
  nni=0;
  while(nni<200){
    //generate a set of neighbors to the current solution
    for(nl_it = neighbor_list->begin(); nl_it != neighbor_list->end(); ++nl_it){
      generateNeighbor(cur_sol, *nl_it);
      calculateObj( *nl_it );
    }
    //sort by objective function
    neighbor_list->sort( objCmp );
    for(nl_it = neighbor_list->begin(); nl_it != neighbor_list->end(); ++nl_it){
      isTabu = 0;
      for(tl_it = tabu_list->begin(); tl_it != tabu_list->end(); ++tl_it){
        if( cmpTabu( *nl_it, *tl_it ) ){
          if( (*nl_it)->obj < best_sol->obj) isTabu = 0;
          else isTabu = 1;
          break;
        }
      }//end for tl_it
      if(isTabu) continue;  // try next neighbor
      cur_sol = *nl_it; //found a good one set it and go on
      break;
    }//end for nl_it
    if(isTabu){ //the back up plan if all neighbors are tabu
      cur_sol = neighbor_list->front();
      std::cout << "warning no non-tabu solution if this happens a lot it could be a problem" << std::endl;
    }
    //current solution to the tabu list.
    tmpSolPtr = tabu_list->back();
    tabu_list->pop_back();
    tabu_list->push_front(tmpSolPtr);
    assignSol(cur_sol, tmpSolPtr);

    if(cur_sol->obj < best_sol->obj){
      nni = 0;
      assignSol(cur_sol, best_sol);
    }
    else{
      ++nni;
    }
    std::cout << "non-improving iterations: " << nni << " best objective: " << best_sol->obj << " current objective: " << cur_sol->obj << std::endl;
  }// end while nni <= nni_limit

  std::cout << "--------------TABU SEARCH RESULT-----------------" << std::endl;
  std::cout << "Initial objective function: " << init_sol->obj << std::endl;
  std::cout << "Best objective function: " << best_sol->obj << std::endl;
  std::cout << std::endl;
  printSol(best_sol);
}

bool pd2_tabuSearch::objCmp(solution_t *A, solution_t *B){
  if(A->obj > B->obj) return 0;
  else if(A->obj < B->obj) return 1;
  else return 0;
}

void pd2_tabuSearch::setTabuList(std::list<solution_t*> *tl){
  tabu_list = tl;
}

void pd2_tabuSearch::setNeighborList(std::list<solution_t*> *nl){
  neighbor_list = nl;
}

void pd2_tabuSearch::setNNILimit(unsigned int nl){
  nni_limit = nl;
}

void pd2_tabuSearch::setBestSol(solution_t *bs){
  best_sol = bs;
}

void pd2_tabuSearch::setInitSol(solution_t *is){
  init_sol = is;
}

void pd2_tabuSearch::setCurrSol(solution_t *cs){
  cur_sol = cs;
}

void pd2_tabuSearch::setTabuCmp( bool (*f)(solution_t *solA, solution_t *solB) ){
  cmpTabu = f;
}

void pd2_tabuSearch::setPrintSol( void (*f)(solution_t *sol) ){
  printSol = f;
}

void pd2_tabuSearch::setCalcObj( void (*f)(solution_t *sol) ){
  calculateObj = f;
}

void pd2_tabuSearch::setGenNeigh( void (*f)(solution_t *sol, solution_t *newSol)  ){
  generateNeighbor = f;
}

void pd2_tabuSearch::setAssign( void (*f)(solution_t *from_sol, solution_t *to_sol)  ){
  assignSol = f;
}

int pd2_tabuSearch::randomIntUniform(int a, int b){
  return(  rand()/(RAND_MAX/(b-a+1)) + a );
}

double pd2_tabuSearch::randomDoubleUniform(double a, double b ){
  return drand48()*(b-a) + a;
}

double pd2_tabuSearch::randomNorm(){
  return(  cos(2*M_PI*(1-drand48()))*sqrt(-2*log(1-drand48()))  );
}

tabusearch.hpp

#ifndef PD2_TABUSEARCH_HPP
#define PD2_TABUSEARCH_HPP

#include<list>
#include<algorithm>
#include<iostream>

class pd2_tabuSearch{
public:
  struct solution_t{
    void *sol;
    double obj;
  };
private:
  std::list<solution_t*> *neighbor_list;
  std::list<solution_t*> *tabu_list;
  solution_t *cur_sol;
  solution_t *init_sol;
  solution_t *best_sol;
  unsigned int nni;
  unsigned int nni_limit;
  void (*printSol)(solution_t *sol); //print some information about a solution
  void (*calculateObj)(solution_t *sol); //calculate an objective function and store in sol
  bool (*cmpTabu)(solution_t *solA, solution_t *solB); //returns 1 if solutions match (for tabu)
  void (*generateNeighbor)(solution_t *sol, solution_t *newSol); //make neighbor list
  void (*assignSol)(solution_t *from_sol, solution_t *to_sol);
  static bool objCmp(solution_t *A, solution_t *B);
public:
  pd2_tabuSearch();
  double initialObjective();
  void search();
  void setTabuList(std::list<solution_t*> *tl);
  void setNeighborList(std::list<solution_t*> *nl);
  void setNNILimit(unsigned int);
  void setBestSol(solution_t *bs);
  void setInitSol(solution_t *is);
  void setCurrSol(solution_t *cs);
  void setTabuCmp( bool (*f)(solution_t *solA, solution_t *solB) );
  void setCalcObj( void (*f)(solution_t *sol) );
  void setGenNeigh( void (*f)(solution_t *sol, solution_t *newSol) );
  void setAssign(void (*f)(solution_t*, solution_t*) );
  void setPrintSol( void (*f)(solution_t *sol) );
  static int randomIntUniform(int a, int b);
  static double randomDoubleUniform(double a, double b);
  static double randomNorm();
};

#endif
