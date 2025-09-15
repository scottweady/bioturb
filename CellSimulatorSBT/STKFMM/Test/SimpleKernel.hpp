#ifndef SIMPLEKERNEL_HPP
#define SIMPLEKERNEL_HPP

//                         3           3           4             4/16/9/7
void StokesSLPVel(double *s, double *t, double *f, double *pvel);
void StokesSLPVelGrad(double *s, double *t, double *f, double *pvelGrad);
void StokesSLTraction(double *s, double *t, double *f, double *traction);
void StokesSLPVelLaplacian(double *s, double *t, double *f, double *pvelLaplacian);

//                         3           3           9             4/16/9/7
void StokesDLPVel(double *s, double *t, double *db, double *pvel);
void StokesDLPVelGrad(double *s, double *t, double *db, double *pvelGrad);
void StokesDLTraction(double *s, double *t, double *db, double *traction);
void StokesDLPVelLaplacian(double *s, double *t, double *db, double *pvelLaplacian);

//                         3           3           1             4/10
void LaplaceSLPGrad(double *s, double *t, double *q, double *pgrad);
void LaplaceSLPGradGrad(double *s, double *t, double *q, double *pgradgrad);

//                         3           3           3             4/10
void LaplaceDLPGrad(double *s, double *t, double *db, double *pgrad);
void LaplaceDLPGradGrad(double *s, double *t, double *db, double *pgradgrad);

//                         3           3           9             10
void LaplaceQPGradGrad(double *s, double *t, double *q, double *pgradgrad);

//
void StokesRegSLVel(double *s, double *t, double *f, double *vel);
void StokesRegSLVelOmega(double *s, double *t, double *f, double *velomega);

//
void StokesRegDLVel(double *s, double *t, double *f, double *vel);
void StokesRegDLVelOmega(double *s, double *t, double *f, double *velomega);

//
void StokesSLRPY(double *s, double *t, double *f, double *vlapv);
void StokesDLRPY(double *s, double *t, double *f, double *vlapv);

void StokesSL(double *s, double *t, double *f, double *v);
void StokesDL(double *s, double *t, double *f, double *v);

//

#endif
