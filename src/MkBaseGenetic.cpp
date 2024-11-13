//---------------------------------------------------------------------------
// This module is core class for Genetic algorithm. The basis of his module
// was built in 1995 to 1996 at Mineral & Petroleum Eng. Dept. of Hanyang
// University. The main work was done by Dr. I.H. Yoo and M.K. Song.
// I.H. Yoo and M.K. Song. habe the copyright of this module, and none can use
// it without notification to the authors.

// It is part of OPTIFACE program which is to be designed to get
// the best parameter using Genetic alogrithm which is one of optimization
// algorithm.
// The work of this project was started by T.H. Ahn, and handed to M.K. Song.
// Copyright (c) 1999 ESCO Consultant Co., Ltd.

#ifndef max
#define max(a, b)  (((a) > (b)) ? (a) : (b))
#endif

#include "MkBaseGenetic.h"

//#include "ga_eval.cpp"

// ----------------------------------------------------------
// INITIALIZE reads in the upper and lower bounds of each
// variable from a text file. Then it randomly generates
// values between these bounds for each gene of each genotype
// in the population. Be aware that the random number
// generator is not seeded using the system time.
//
// The format of the input file "gadata.txt"is as follows:
// Var1_lower_bound  Var1_upper-bound
// Var2_lower_bound  Var2_upper-bound
// etc.
// ----------------------------------------------------------
int MAXGENS;		// Maximum number of generations
int Generation;     // Current generation number

MkBaseGenetic::MkBaseGenetic()
{
     MAXGENS = 0;
     NumOfVar = 0;
     MaxCrom = 0;
     NumOfPop = 0;
     IsMax = 0;
     Pc = 0;
     Pm = 0;

     Gen  = NULL;
     Upper = NULL;
     Lower = NULL;
     Prob = NULL;
     ProbSum = NULL;
     Fitness = NULL;
}

MkBaseGenetic::~MkBaseGenetic()
{
     if (Gen) delete[] Gen;
     if (Upper) delete[] Upper;
     if (Lower) delete[] Lower;
     if (Prob) delete[] Prob;
     if (ProbSum) delete[] ProbSum;
     if (Fitness) delete[] Fitness;
}

bool MkBaseGenetic::Initialize(int num_of_var, int num_of_pop, bool is_max, float pc,
     float pm, int max_gen)       // memory alloc
{

     int i,j;
     int pop_size,gen_max,numvar;

     pop_size = num_of_pop;
     gen_max = max_gen;
     numvar = num_of_var;

     if(pop_size <= 0) {
       //ShowMessage("Pop size is zero or wrong data file for GAWIN\n");
        return false;
     }

     if(gen_max <= 0) {
       //ShowMessage("Max Generation is zero or wrong data file for GAWIN\n");
        return false;
     }

     if(numvar <= 0 || numvar > 5 ) {
       //ShowMessage("Numof of variable is zero or wrong data file for GAWIN\n"
       //       "or Exceed the Maximum variable number(=5) for evaluation copy\n");
        return false;
     }

     MAXGENS = gen_max;

     NumOfVar = numvar;
     MaxCrom = NumOfVar*20;
     NumOfPop = pop_size;
     IsMax = is_max;
     Pc = pc;
     Pm = pm;

     Gen  = new char[NumOfPop][100];
     if(Gen  == NULL ) return false;

     Upper = new double[NumOfVar];
     if(Upper == NULL ) return false;

     Lower = new double[NumOfVar];
     if(Lower == NULL ) return false;

     Prob = new double[NumOfPop];
     if(Prob == NULL ) return false;

     ProbSum = new double[NumOfPop];
     if(ProbSum == NULL ) return false;

     Fitness = new double[NumOfPop];
     if(Fitness == NULL ) return false;

     // Using the lower and upper bounds, randomly assign a value to
     // each gene of each genotype.
}

void MkBaseGenetic::Initialize(int num, float low, float high)
{
     Lower[num] = low;
     Upper[num] = high;
}

// ----------------------------------------------------------
// RandVal returns a pseudorandom value between the supplied
// low and high arguments.
// ----------------------------------------------------------

double MkBaseGenetic::RandVal(double low, double high)
{
	double val;

	val = ((double)(rand() % 10000) / 10000.0) * (high - low) + low;
	return (val);
}

// ----------------------------------------------------------
// Evaluate is a user defined function. For each problem you
// solve using GA, you will modify this function and
// recompile the code. The simple function "a^2 - ab + c" is
// supplied here for illustration purposes.
// ----------------------------------------------------------

void MkBaseGenetic::Evaluate()
{
	int    mem;
	int    i;

//	Best = 0;

	for (mem = 0; mem < NumOfPop; mem++)
            Fitness[mem] = Eval(mem);

	if ( IsMax ) {
		for (mem = 0; mem < NumOfPop; mem++)
			if( BestVal < Fitness[mem] ) {
				BestVal = Fitness[mem];
                                BestPop = mem;
				for ( int i = 0 ; i < MaxCrom ; i++ )
					BestGen[i] = Gen[mem][i];
			}
	}
	else {
		for (mem = 0; mem < NumOfPop; mem++)
			if( BestVal > Fitness[mem] ) {
				BestVal = Fitness[mem];
                                 BestPop = mem;
				for ( int i = 0 ; i < MaxCrom ; i++ )
					BestGen[i] = Gen[mem][i];
			}
	}
	double worst_val = Fitness[0];
	int W;W=0;

	for (i = 0; i < NumOfPop; i++)
	{
		if (IsMax) {
			if (worst_val > Fitness[i]) {
				worst_val =  Fitness[i];
				W = i;
			}
		}
		else {
			if (worst_val < Fitness[i]) {
				worst_val =  Fitness[i];
				W = i;
			}
		}
	}

	if(W>0&&W<NumOfPop) {
		for ( int k = 0 ; k < MaxCrom ; k++) {
			Gen[W][k] = BestGen[k];
		}
		Fitness[W] = BestVal;

	}
}

// ----------------------------------------------------------
// This is an implementation of standard single point
// crossover. Many other crossover operators developed
// specifically for real-coded GA's may give better results.
// For simplicity, only the single point crossover is
// shown here.
// ----------------------------------------------------------

double MkBaseGenetic::Decode(int pop,int var_num)
{
	long sunsea=0;
	double value;
	long Total = (long)(pow((double)2,(double)20)-1);

	for ( int i = 0+var_num*20 ; i < 20+var_num*20 ; i ++ )
		sunsea+= (long)pow(double(2),double(var_num*20+19-i))*Gen[pop][i];

	value = Lower[var_num] + (Upper[var_num]-Lower[var_num])*sunsea/Total;
	return value;
}

double MkBaseGenetic::Decode(char gen[100],int var_num)
{
	long sunsea=0;
	double value;
	long Total = (long)(pow((double)2,(double)20)-1);

	for ( int i = 0+var_num*20 ; i < 20+var_num*20 ; i ++ )
		sunsea+= (long)pow(double(2),double(var_num*20+19-i))*gen[i];

	value = Lower[var_num] + (Upper[var_num]-Lower[var_num])*sunsea/Total;
	return value;
}

void MkBaseGenetic::SetWheel(void)
{
	double sum;
	sum = 0;
	double min;
        int  i;
	min = Fitness[0];

	for ( i = 0 ; i < NumOfPop ; i++ )
		min = min > Fitness[i] ? Fitness[i] : min;
	if( min < 0 ) {

		for ( i = 0 ; i < NumOfPop ; i++ ) {
			sum += Fitness[i]-2*min;
		}
		double mean = sum / NumOfPop;

		for ( i = 0 ; i < NumOfPop ; i++ ) {
			if( IsMax ) Prob[i] = (Fitness[i]-2*min)/sum;
			else Prob[i] = (2*mean-Fitness[i]+2*min)/sum;
		}
	}
	else {
		for ( i = 0 ; i < NumOfPop ; i++ )
			sum += Fitness[i];

		double mean = sum / NumOfPop;

		for ( i = 0 ; i < NumOfPop ; i++ ) {
			if( IsMax ) Prob[i] = Fitness[i]/sum;
			else Prob[i] = (2*mean-Fitness[i])/sum;
		}
	}

	ProbSum[0] = Prob[0];

	for ( i = 1 ; i < NumOfPop ; i++ )
		ProbSum[i] = Prob[i]+ProbSum[i-1];
}

void MkBaseGenetic::Roulette()
{
	double *R;
        int  i,j,k;
	R = new double[NumOfPop];

	if(R==NULL) return;

	for ( i = 0 ; i < NumOfPop ; i++ )
		R[i] = RandVal(0,1);

	for (i = 0 ; i < NumOfPop ; i++ ) {
		for (j = 0 ; j < NumOfPop ; j++ )
			if( R[i] < ProbSum[j]) break;

		for ( k = 0 ; k < MaxCrom ; k++) {
			Gen[i][k] = Gen[j][k];
		}

	}

	delete[] R;
}

void MkBaseGenetic::CrossOver(void)
{
	float *p;
	int *q;
	int i,j;
	char tmp;
	p = new float[NumOfPop];
	q = new int[NumOfPop];

	if(p==NULL||q==NULL) return;

	for(i = 0 ; i < NumOfPop;i++ )
		p[i] = RandVal(0,1);

	for(i = 0,j=0 ; i < NumOfPop;i++ )
		if(p[i]<Pc) {q[j] = i;j++;}

	int NumOfCross = j;

	int I,J;
	I = (int)RandVal(5,15);
	J = (int)RandVal(25,35);

	for ( i = 0 ; i < NumOfCross/2 ;i+=2) {
		for ( j = I ; j < J ; j++ ) {
			tmp = Gen[q[i]][j];
			Gen[q[i]][j] = Gen[q[i+1]][j];
			Gen[q[i+1]][j] = tmp;
		}
	}
	delete[] q;
	delete[] p;
}

// ----------------------------------------------------------
// This is an implementation random mutation. A value
// selected for mutation is replaced by a randomly generated
// number between that variable's lower and upper bounds.
// As is the case with crossover, other mutations developed
// for real-coded GA's can be added.
// ----------------------------------------------------------

void MkBaseGenetic::Mutate(void)
{
	int *p;
	int i,j;
	int NumBit;
	NumBit = (int)(NumOfPop*MaxCrom*Pm*3);
	p = new int[NumBit];

	if(p==NULL) return;

	for(i = 0,j = 0 ; i < NumOfPop*MaxCrom && j < NumBit;i++ )
		if (RandVal(0,1)<Pm) {p[j] = i;j++;}

	int NumOfMutate = j;

	for ( i = 0 ; i < NumOfMutate ;i++) {
		int Q,R;
		Q = p[i]/MaxCrom;
		R = p[i]%MaxCrom;
		Gen[Q][R] = Gen[Q][R] ? 0 : 1;
	}
	delete[] p;
}

//	----------------------------------------------------------
//	Report progress of the simulation. Data is dumped to a log
//	file in comma seperated value format which can be imported
//	and graphed using any commericial spreadsheet pcckage.
//	----------------------------------------------------------
/*
void MkBaseGenetic::Report(TObject *Sender)
{
	int i,I=0,W=0;
	double best_val;     // Best population fitness
	static double old_best;
	double worst_val;     // Best population fitness
	double X[5];
        char str[256],str2[256];

	// Calculate the summed population fitness and the sum of the
	// squared individual fitness

	worst_val = best_val = Fitness[0];
	W = 0;

	for (i = 0; i < NumOfPop; i++) {
            if (IsMax) {
               if (worst_val > Fitness[i]) {
                  worst_val =  Fitness[i];
                  W = i;
               }
            }
            else {
               if (worst_val < Fitness[i]) {
                  worst_val =  Fitness[i];
                  W = i;
               }
            }
	}

	if( W >= 0 && W < NumOfPop ) {
            for ( int k = 0 ; k < MaxCrom ; k++) {
                Gen[W][k] = BestGen[k];
            }
            Fitness[W] = BestVal;
	}
	if( W < 0 || W >= NumOfPop ) {
            for ( int k = 0 ; k < MaxCrom ; k++) {
                Gen[0][k] = BestGen[k];
            }

            Fitness[0] = BestVal;
	}

	for (i = 0; i < NumOfPop; i++)
	{
            if ( IsMax ) {
               if (best_val < Fitness[i]) {
                  best_val =  Fitness[i];
                  I = i;
               }
            }
            else {
               if (best_val > Fitness[i]) {
                  best_val =  Fitness[i];
                  I = i;
               }
            }
	}

	for ( i = 0 ; i < NumOfVar ; i ++ )
            X[i] = Decode(I,i);

        if( best_val != old_best ) {
            sprintf(str,"\n%5d> %10.5f | ",Generation,best_val);
            for ( i = 0 ; i < NumOfVar ; i++ ) {
                sprintf(str2," %10.5f ",X[i]);
                strcat(str,str2);
            }
            old_best = best_val;
        }
        if (String(Sender->ClassName()) == String("TMemo")) {
           TMemo *Memo1=(TMemo*)Sender;
           Memo1->Lines->Add(str);
        }
}
*/

float  MkBaseGenetic::GetFitness(int i)
{
    if (i>=0&&i<NumOfPop) return Fitness[i];
    else if (IsMax) return -100000;
    else if (!IsMax) return 100000;
};


double MkBaseGenetic::Eval(int pop)
{
    return Eval(Gen[pop]);
}

double MkBaseGenetic::Eval(char gen[100])
{
    return 0;
}


float MkBaseGenetic::frandom(float f)
{
  //    return random((int)f*1000)/1000.0;
  return f*float(rand())/RAND_MAX;
}

//---------------------------------------------------------------------------
#pragma package(smart_init)
