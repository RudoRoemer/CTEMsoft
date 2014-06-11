#include "data.h" 
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "BridgeTest.h"


static Data* g_Data = NULL;


// Callback to allocate Memory
void CppProg(int FZcnt, double* &outp)
{

  outp = static_cast<double*>(malloc(FZcnt * 3 * sizeof(double)));
  ::memset(outp, 0xAB, FZcnt * 3 * sizeof(double));
  printf("outp: %p \n", outp);
  g_Data->setDataPointer(outp);
  g_Data->setFZCount(FZcnt);
}


Data::Data() :
m_PointGroup(0),
m_NumSteps(0),
m_Data(NULL),
m_FZCount(0)
{
	g_Data = this;
}

Data::~Data()
{
	if(NULL != m_Data) {  free(m_Data);  }
}

void Data::setFZCount(int fzcount)
{
	m_FZCount = fzcount;
}

int Data::getFZCount()
{
	return m_FZCount;
}

void Data::setDataPointer(double* d)
{
	m_Data = d;
}
double* Data::getDataPointer()
{
	return m_Data;
}


void Data::setPointGroup(int pgnum)
{
	m_PointGroup = pgnum;
}
int Data::getPointGroup()
{
	return m_PointGroup;
}


void Data::setNumSteps(int nsteps)
{
 m_NumSteps = nsteps;
}
int Data::getNumSteps()
{
return m_NumSteps;
}

void Data::execute()
{
 	TheCStruct_t mycstruct; 
	mycstruct.i = 11;
	mycstruct.j = 12;
	mycstruct.s = 2;
	mycstruct.m = 1.5;
	mycstruct.d = 2.5;

	sampler(&m_PointGroup, &m_NumSteps, &mycstruct, &CppProg);

}