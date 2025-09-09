#ifndef PARAMETER
#define PARAMETER
struct parameter{
	double t_oxide; //oxide thickness
	double e_oxide;
	double e_si;
	double Na;
	double mobility;
	double V;
	double size;
};
extern struct parameter mos;
#endif
