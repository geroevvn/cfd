#include "CellFluidDynamicsProps.h"

CellFluidDynamicsProps::CellFluidDynamicsProps(double ro, double u, double v, double w, double E,
												double P, double gamma) {
	this->ro = ro;
	this->ru = ro * u;
	this->rv = ro * v;
	this->rw = ro * w;
	this->rE = ro * E;
	this->P = P;
	this->gamma = gamma;
}

CellFluidDynamicsProps::CellFluidDynamicsProps() {
}


CellFluidDynamicsProps::~CellFluidDynamicsProps() {
}

double CellFluidDynamicsProps::calc_rE(double ro, double P, double u, double v, double w, double gamma)
{
	return P / (gamma - 1) + ro * 0.5 *(u*u + v*v + w*w);
}

double CellFluidDynamicsProps::calc_P(double ro, double rE, double ru, double rv, double rw, double gamma)
{
	return (rE - 0.5 *(ru*ru + rv*rv + rw*rw) / ro) * (gamma - 1);
}

double CellFluidDynamicsProps::ro_Roe_Avg(double ro1, double ro2)
{
	return sqrt(ro1*ro2);
}

double CellFluidDynamicsProps::vel_Roe_Avg(double ro1, double ro2, double u1, double u2)
{
	double ro1_sqrt = sqrt(ro1);
	double ro2_sqrt = sqrt(ro2);

	return ( ro1_sqrt*u1 + ro2_sqrt*u2 ) / ( ro1_sqrt + ro2_sqrt );
}

double CellFluidDynamicsProps::H_Roe_Avg(double ro1, double ro2, double H1, double H2)
{
	double ro1_sqrt = sqrt(ro1);
	double ro2_sqrt = sqrt(ro2);

	return (ro1_sqrt*H1 + ro2_sqrt*H2) / (ro1_sqrt + ro2_sqrt);
}

double CellFluidDynamicsProps::rE_Roe_Avg(double ro, double u, double v, double w, double gamma, double H)
{
	return ro * (H + (gamma - 1)*(u*u + v*v + w*w)/2) / gamma;
}

double CellFluidDynamicsProps::P_Roe_Avg(double ro, double u, double v, double w, double gamma, double H)
{
	return ro * (1 - 1.0/gamma) * (H - (u*u + v*v + w*w)/2);
}

double CellFluidDynamicsProps::calc_Roe_Avg(double& u, double& v, double& w, double& H, double& c, double& GAMMA, const CellFluidDynamicsProps& cfdp1, const CellFluidDynamicsProps& cfdp2)
{
	double P;
	double ro;
	double E;

	double u1 = cfdp1.ru / cfdp1.ro;
	double v1 = cfdp1.rv / cfdp1.ro;
	double w1 = cfdp1.rw / cfdp1.ro;

	double u2 = cfdp2.ru / cfdp2.ro;
	double v2 = cfdp2.rv / cfdp2.ro;
	double w2 = cfdp2.rw / cfdp2.ro;

	GAMMA = cfdp1.gamma;
	double AGAM = GAMMA - 1.0;

	double fSB = sqrt( cfdp1.ro );
	double fSE = sqrt( cfdp2.ro );
	double fS_ = 1.0 / ( fSB + fSE );

	ro = fSB * fSE;

	u = ( fSB * u1 + fSE * u2 ) * fS_;
	v = ( fSB * v1 + fSE * v2 ) * fS_;
	w = ( fSB * w1 + fSE * w2 ) * fS_;

	double EB = cfdp1.P / (cfdp1.ro * AGAM);
	double EE = cfdp2.P / (cfdp2.ro * AGAM);

	double HB = EB + 0.5 * (u1*u1 + v1*v1 + w1*w1) + cfdp1.P / cfdp1.ro;
	double HE = EE + 0.5 * (u2*u2 + v2*v2 + w2*w2) + cfdp2.P / cfdp2.ro;

	double HI = ( fSB * HB + fSE * HE ) * fS_;

	P = (HI - 0.5 * (u*u + v*v + w*w) ) * ro * AGAM / GAMMA;
	E = P / (ro * AGAM) + 0.5 * (u*u + v*v + w*w);

	c = sqrt(GAMMA * P / ro);
	H = E + P / ro;
}
