// ***********************************
// This program is developed for the reserch of freezing and melting of ice
// This is basic test bed class for upcomming SPH simulator
// ***********************************

//---------------------------------------------------------------------------
#ifndef SPHH
#define SPHH
//---------------------------------------------------------------------------

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <GL/glut.h>
#include "rvgs.h"
#include "MkMatrix.h"
#include "MkSmoothFunc.h"
#include "MkParticle.h"
#include "MkGrid.h"

class MkSPH
{
protected:
	int NumOfParticles;
	int GI, GJ;
	float LowX, LowY, UpX, UpY;
	float RefDensity, Gamma_P, SoundSpeed, XGravity, YGravity;
	float Time, Iteration, DeltaTime;
	float SmoothLen;
	MkParticles Particles;
	MkSmoothFunc SmoothFunc, SmoothFuncTemp;
	MkSPHGrids SPHGrids;
	MkInt Grids;

public:
	MkSPH();
	~MkSPH();

	bool Initialize();											 // memory allocation
	bool Initialize(int np, int gi, int gj); // memory allocation
	void Clear();
	void SetNumOfParticles(int n) { NumOfParticles = n; }
	void SetGIGJ(float gi, float gj)
	{
		GI = gi;
		GJ = gj;
	}
	void SetGI(float gi) { GI = gi; }
	void SetGJ(float gj) { GJ = gj; }
	void SetBound(float lx, float ly, float ux, float uy)
	{
		LowX = lx;
		LowY = ly;
		UpX = ux;
		UpY = uy;
	}
	void SetRefDensity(float refden) { RefDensity = refden; }
	void SetGammaPolytropic(float gp) { Gamma_P = gp; }
	void SetSoundSpeed(float ss) { SoundSpeed = ss; }
	void SetSmoothLen(float sm)
	{
		SmoothLen = sm;
		SmoothFunc.SetupSmoothFunc(smtGaussian, 2, SmoothLen, 1);
		SmoothFuncTemp.SetupSmoothFunc(smtGaussian, 2, SmoothLen, 1);
	}
	void SetGravity(float xg, float yg)
	{
		XGravity = xg;
		YGravity = yg;
	}
	void SetDeltaTime(float dt) { DeltaTime = dt; }

	bool GenParticle();
	bool GenParticle2();
	MkParticles &GetParticles() { return Particles; }
	int GetNumOfParticles(void) { return NumOfParticles; }
	float GetGI(void) { return GI; }
	float GetGJ(void) { return GJ; }
	float GetRefDensity() { return RefDensity; }
	float GetGammaPolytropic() { return Gamma_P; }
	float GetSoundSpeed() { return SoundSpeed; }
	float GetSmoothLen() { return SmoothLen; }
	void GetGravity(float &xg, float &yg)
	{
		xg = XGravity;
		yg = YGravity;
	}
	float GetDeltaTime(void) { return DeltaTime; }
	float GetTime(void) { return Time; }

	bool InitSPHGrid();
	bool InitAccel();
	bool InitVelocity();
	bool InitPosition();
	bool UpdateSPHGrid();
	bool ComputeDensity();
	bool ComputeDensityRate();
	bool ComputePressure();
	bool ComputeInternalForce2();
	bool ComputeInternalForce();
	bool ComputeViscousForce2();
	bool ComputeArtificialViscousForce2();
	bool ComputeViscousForce();
	bool ComputeVelocity();
	bool ComputeAccel();
	bool ComputePosition();
	bool UpdatePosition();

public:
	bool Run();
};

#endif
