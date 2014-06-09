#include <iostream>
#include <string>
#include <sstream>
#include "BridgeTest.h"


double* g_out = NULL;

void CppProg(int FZcnt, double* &outp)
{

  outp = static_cast<double*>(malloc(FZcnt * 3 * sizeof(double)));
  ::memset(outp, 0xAB, FZcnt * 3 * sizeof(double));
  printf("outp: %p \n", outp);
  g_out = outp;
}
//
// void StrProg(const char* str)
// {
//
//     for(int i = 0; i < 40; i++)
//     {
//         std::cout << "I=" << i << "   " << (int)(str[i]) << std::endl;
//     }
// 	std::cout << str << std::endl;
//
// }


int main(int argc, char* argv[])
{
	std::cout << "Fortran to C++ Bridge..." << std::endl;
	std::string mystr;
	int pgnum=0;
 	int nsteps=0;
 	double* outp = NULL;

  	std::cout << "Enter the point group number :  ";
  	std::cin >> pgnum;


  	std::cout << "Enter the number of intervals along the cube semi-edge length : ";
  	std::cin >> nsteps;


 	sampler_(&pgnum, &nsteps, &CppProg);


 	std::cout << "Done" << std::endl;
 	std::cout << *(g_out+0) << std::endl;
	std::cout << *(g_out+1) << std::endl;
	std::cout << *(g_out+8069) << std::endl;

 	//std::cout << "Test" << std::endl;
//  	nmltypefile nmlfile;
//
//  	nmlfile.sig=70.0;




 	//std::string inputFile = "test.nml";
//  	domcsimulation_(&nmlfile);
 	//teststr_(&StrProg);

	//std::cout << result << std::endl;

	return EXIT_SUCCESS;
}

