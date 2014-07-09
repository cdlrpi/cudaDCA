//	classes.h
//
//This file contains the body and joint classes used to describe the multibody system.

//InitBody:
//	Body class used only during initialization
//	The bodies are separated in this way to avoid carrying unnecessary 
//	information throughout the DCA algorithm.  The mass, inertia, and length
//	of each body is not relevant once the initial zetas are found
class InitBody
{
	public:
		double m; //mass
		double I[3][3]; //inertia tensor
		double l; //length
};

//Body:
//	Body class used throughout the DCA algorithm.  This class has no need for m, I and l
class Body
{
	public:

		//Zeta values
		double z11[6][6];
		double z12[6][6];
		double z13[6];
		double z21[6][6];
		double z22[6][6];
		double z23[6];

		//Xinv is an intermediate quantity that is stored because of its
		//presence in assembly and disassembly
		double Xinv[5][5];

};

//Joint:
//	Joint class that initiates a 2d pendulum joint
//	NOTE:  This class is not used in this version
class Joint
{
	public:
		Joint(); //Initialization Function prototype
};
