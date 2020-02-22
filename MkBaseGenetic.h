// ***********************************
// This program is developed for the reserch of Dr. Yoo, by both him and me.
// At his thesis, this algorithm was used to find optimum parameter. It was
// similar to least square method.
// This class was first used to develope GAWIN 1.0 at 1996.3.
// Also, it is modified for C++ Builder at 99.2.26 - 99.6
// This code is used for optimizing the tunnel face, to minimize the
// excavation cost, in the field of tunnel engineering.
// ***********************************

//---------------------------------------------------------------------------
#ifndef BaseGeneticH
#define BaseGeneticH
//---------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include <ctype.h>

// Change any of these parameters to match your needs

#define MAXTRY 10
extern int MAXGENS;		// Maximum number of generations
extern int Generation;     // Current generation number

// Each genotype is a member of the population

class MkBaseGenetic {
protected:
	int NumOfVar;
	int NumOfPop;
	bool IsMax;		   // Minimization = false, Maximize = true
	float Pc ;		   // Percent of crossover
	float Pm ;		   // Percent of mutation
	double *Upper;             // The genotype's upper bound
	double *Lower;             // The genotype's lower bound
	double *Prob;              // What's this?
	double *ProbSum;           // What's this too?
	char   (*Gen)[100];          // A string of bits makes a genetype

	double *Fitness;           // The genotype's fitness

	double BestVal;
        int    BestPop;
	char   BestGen[100];

	int    MaxCrom;

public:
        MkBaseGenetic();
	~MkBaseGenetic();

        virtual bool Initialize( int num_of_var,int num_of_pop, bool is_max,
                         float pc, float pm, int max_gen); // memory allocation
        virtual void Initialize( int var_num, float low, float high); // boundary of var

	virtual double Eval(int pop);  // pure virtual routine
	virtual double Eval(char gen[100]);  // pure virtual routine
	virtual double Decode(int pop,int var_num);  // pure virtual routine
	virtual double Decode(char gen[100],int var_num);  // pure virtual routine

	double RandVal(double low, double high);
	void Evaluate();
	void CrossOver(void);
	void Mutate(void);
	//	void Report(TObject *Sender);
	void SetWheel(void);
	void Roulette(void);
	int GetNumOfVar(){return NumOfVar;};
	double GetBestVal(){return BestVal;};
        char * GetBestGen(){return BestGen;};
        int    GetBestPop(){return BestPop;};
        float  GetFitness(int i);
        int    GetNumOfPop(){return NumOfPop;};
        float frandom(float f);
public:
};

#endif
