/*                                   classes.h

	This file contains the body and joint classes used to describe the multibody system.
*/
#ifndef __CLASSES_H_INCLUDED__
#define __CLASSES_H_INCLUDED__
//Body class used during initialization
class InitBody
{
	public:
		float m; //mass
		float I[3][3]; //inertia tensor
		float l; //length
};

//Body class used only during assembly and disassembly.  This class has no need for m, I and l
class Body
{
	public:

		//Zeta values
		float z11[6][6];
		float z12[6][6];
		float z13[6];
		float z21[6][6];
		float z22[6][6];
		float z23[6];
		float Xinv[5][5];

};

//Joint class that initiates to a 2d pendulum joint
class Joint
{
	public:
		Joint(); //Initialization Function prototype
};



#endif
