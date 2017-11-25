#include "ShellEnergy.h"
#include "MaterialParameters.h"
#include "Geometry.h"
#include "Eigen/Dense"
#include <iostream>

using namespace Eigen;

scalar triangleStretchingEnergy(
    const Vector3m &p0,
    const Vector3m &p1,
    const Vector3m &p2,
    const Matrix2m &abar,
    scalar h,
    const MaterialParameters &params,
    Eigen::Matrix<scalar, 9, 1> *dEnergy,
    Eigen::Matrix<scalar, 9, 9> *hEnergyExact,
    Eigen::Matrix<scalar, 9, 9> *hEnergyInexact)
{
    Matrix2m a = firstFundamentalForm(p0, p1, p2);
    Matrix2m abarinv = abar.inverse();
    Matrix2m M = abarinv*(a - abar);
    Matrix <scalar, 3, 2 > r;
    r.col(0) = p1 - p0;
    r.col(1) = p2 - p0;

    scalar dA = 0.5 * sqrt(abar.determinant());

    scalar energy = h / 4.0 * (params.LameAlpha() / 2.0 * M.trace()*M.trace() + params.LameBeta() * (M*M).trace()) * dA;

    Vector2m wi[3];
    wi[0] << -1, -1;
    wi[1] << 1, 0;
    wi[2] << 0, 1;

    Eigen::Matrix<scalar, 9, 1> deriv1;
    Eigen::Matrix<scalar, 9, 1> deriv2;

    if (dEnergy || hEnergyInexact)
    {
        
        for (int i = 0; i < 3; i++)
        {
            deriv1.segment<3>(3 * i) = 2.0*r*abarinv*wi[i];
            deriv2.segment<3>(3 * i) = 2.0*r*abarinv*M.transpose()*wi[i];
        }
        
        if(dEnergy)
            *dEnergy = h*params.LameAlpha()* dA / 4.0 * M.trace() * deriv1 + h*params.LameBeta() * dA / 2.0 * deriv2;        

        if (hEnergyInexact)
        {
            double denom = (M*M).trace();
            if (denom == 0.0)
                hEnergyInexact->setZero();
            else
            {
                *hEnergyInexact = h*params.LameAlpha()*dA / 4.0 * deriv1*deriv1.transpose() + h*params.LameBeta()*dA / 2.0 / denom * deriv2*deriv2.transpose();
            }
        }
    }
    

    if (hEnergyExact)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                Matrix3m hE;
                hE.setIdentity();
                scalar hscale = wi[j].transpose() * (params.LameAlpha()*M.trace() * abarinv + 2.0*params.LameBeta()*abarinv * M.transpose()) * wi[i];
                hscale *= h / 2.0 * dA;
                hE *= hscale;

                Matrix3m hE2 = r*abarinv*wi[i] * wi[j].transpose()*abarinv*r.transpose();
                hE2 *= h * (params.LameAlpha() + params.LameBeta()) * dA;

                Matrix3m hE3 = r*abarinv*r.transpose();
                scalar prod = (wi[j].transpose()*abarinv*wi[i]);
                hE3 *=  prod* h*params.LameBeta()*dA;

                hEnergyExact->block<3, 3>(3 * i, 3 * j) = hE + hE2 + hE3;
            }
        }
    }

    return energy;
}

scalar triangleBendingEnergy(
    const Vector3m &p0,
    const Vector3m &p1,
    const Vector3m &p2,
    const Vector3m *q0,
    const Vector3m *q1,
    const Vector3m *q2,
    const Matrix2m &abar,
    const Matrix2m &bbar,
    scalar h,
    const MaterialParameters &params,
    Matrix2m *a,
    Matrix2m *b,
    Eigen::Matrix<scalar, 18, 1> *dEnergy,
    Eigen::Matrix<scalar, 18, 18> *hEnergyInexact)
{
    Vector3m mc = (p1 - p0).cross(p2 - p0);
    Vector3m fi[3];
    if (q0)
    {
        fi[0] = (p2 - *q0).cross(p1 - *q0);
    }

    if (q1)
    {
        fi[1] = (p0 - *q1).cross(p2 - *q1);
    }

    if (q2)
    {
        fi[2] = (p1 - *q2).cross(p0 - *q2);
    }

    const Vector3m *qi[3];
    qi[0] = q0;
    qi[1] = q1;
    qi[2] = q2;

    Vector3m mi[3];
    Vector3m ni[3];
    for (int i = 0; i < 3; i++)
    {
        mi[i] = mc / mc.norm();
        if(qi[i])
            mi[i] += fi[i] / fi[i].norm();
        ni[i] = mi[i] / mi[i].norm();
    }

    *a = firstFundamentalForm(p0, p1, p2);
    Matrix2m abarinv = abar.inverse();
    Matrix<scalar, 3, 2> r;
    r.col(0) = p1 - p0;
    r.col(1) = p2 - p0;

    Matrix<scalar, 3, 2> n;
    n.col(0) = ni[0] - ni[1];
    n.col(1) = ni[0] - ni[2];
    
    *b = r.transpose()*n + n.transpose()*r;
    Matrix2m M = abarinv*(*b - bbar);
    scalar dA = 0.5 * sqrt(abar.determinant());

    scalar energy = h*h*h / 12.0 * (params.LameAlpha() / 2.0 * M.trace()*M.trace() + params.LameBeta() * (M*M).trace()) * dA;

    Vector2m wi[3];
    wi[0] << -1, -1;
    wi[1] << 1, 0;
    wi[2] << 0, 1;

    Vector3m pi[3];
    pi[0] = p0;
    pi[1] = p1;
    pi[2] = p2;

    if (dEnergy || hEnergyInexact)
    {
        Matrix<scalar, 18, 1> alphav;
        Matrix<scalar, 18, 1> betav;
        alphav.setZero();
        betav.setZero();
        scalar denom = 0.5 / sqrt((M*M).trace());
        for (int i = 0; i < 3; i++)
        {
            Vector3m alphadEi = 2.0*n*(abarinv*wi[i]);
            Vector3m betadEi = 4.0*denom*n*abarinv*(M.transpose()*wi[i]);
            alphav.segment<3>(3 * i) += alphadEi;
            betav.segment<3>(3 * i) += betadEi;

            Matrix3m Id;
            Id.setIdentity();
            Vector3m alphancommon = -2.0*(1.0 / mi[i].norm()) * (Id - ni[i] * ni[i].transpose()) * r*(abarinv*wi[i]);
            Vector3m betancommon = -4.0*denom * (1.0 / mi[i].norm()) * (Id - ni[i] * ni[i].transpose()) * r*abarinv*(M.transpose()*wi[i]);
            
            Matrix3m mcmat = 1.0 / mc.norm() * (Id - mc*mc.transpose() / mc.squaredNorm());
            Vector3m alphaterm1 = crossMatrix(pi[(i + 2) % 3] - pi[i])*mcmat*alphancommon;
            alphav.segment<3>(3 * ((i + 1) % 3)) += alphaterm1;
            alphav.segment<3>(3 * i) -= alphaterm1;

            Vector3m betaterm1 = crossMatrix(pi[(i + 2) % 3] - pi[i])*mcmat*betancommon;
            betav.segment<3>(3 * ((i + 1) % 3)) += betaterm1;
            betav.segment<3>(3 * i) -= betaterm1;

            Vector3m alphaterm2 = -crossMatrix(pi[(i + 1) % 3] - pi[i])*mcmat*alphancommon;
            alphav.segment<3>(3 * ((i + 2) % 3)) += alphaterm2;
            alphav.segment<3>(3 * i) -= alphaterm2;

            Vector3m betaterm2 = -crossMatrix(pi[(i + 1) % 3] - pi[i])*mcmat*betancommon;
            betav.segment<3>(3 * ((i + 2) % 3)) += betaterm2;
            betav.segment<3>(3 * i) -= betaterm2;

            if (qi[i])
            {
                Matrix3m fmat = 1.0 / fi[i].norm() * (Id - fi[i] * fi[i].transpose() / fi[i].squaredNorm());
                Vector3m alphaterm1 = crossMatrix(pi[(i + 1) % 3] - *qi[i])*fmat*alphancommon;
                alphav.segment<3>(3 * ((i + 2) % 3)) += alphaterm1;
                alphav.segment<3>(9 + 3 * i) -= alphaterm1;
                
                Vector3m betaterm1 = crossMatrix(pi[(i + 1) % 3] - *qi[i])*fmat*betancommon;
                betav.segment<3>(3 * ((i + 2) % 3)) += betaterm1;
                betav.segment<3>(9 + 3 * i) -= betaterm1;

                Vector3m alphaterm2 = -crossMatrix(pi[(i + 2) % 3] - *qi[i])*fmat*alphancommon;
                alphav.segment<3>(3 * ((i + 1) % 3)) += alphaterm2;
                alphav.segment<3>(9 + 3 * i) -= alphaterm2;

                Vector3m betaterm2 = -crossMatrix(pi[(i + 2) % 3] - *qi[i])*fmat*betancommon;
                betav.segment<3>(3 * ((i + 1) % 3)) += betaterm2;
                betav.segment<3>(9 + 3 * i) -= betaterm2;
            }
        }

        if ((M*M).trace() == 0.0)
            betav.setZero();

        if (dEnergy)
        {
            dEnergy->setZero();
            *dEnergy = h*h*h / 12.0*dA*(params.LameAlpha()*M.trace()*alphav + 2.0*params.LameBeta()*sqrt((M*M).trace())*betav);
        }
        if (hEnergyInexact)
        {
            hEnergyInexact->setZero();
            *hEnergyInexact = h*h*h / 12.0 * dA*(params.LameAlpha()*alphav*alphav.transpose() + 2.0*params.LameBeta()*betav*betav.transpose());
        }
    }
    
    return energy;
}

scalar shellEnergy(
    const MatrixX3m &V,
    const Eigen::MatrixX3i &F,
    const Eigen::MatrixX3i &faceWings,
    const std::vector<Matrix2m> &abars,
    const std::vector<Matrix2m> &bbars,
    const VectorXm &faceThicknesses,
    const MaterialParameters &params,
    std::vector<Matrix2m> *acurrent,
    std::vector<Matrix2m> *bcurrent,
    VectorXm *dEnergy,
    std::vector<Eigen::Triplet<scalar> > *hEnergyExact,
    std::vector<Eigen::Triplet<scalar> > *hEnergyInexact,
    VectorXm *triangleEnergies)
{
    scalar result=0;
    int nverts = (int)V.rows();
    int nfaces = (int)F.rows();

    if (dEnergy && dEnergy->size() != 3*nverts)
    {
        dEnergy->resize(nverts * 3);
        dEnergy->setZero();
    }

    if (triangleEnergies && triangleEnergies->size() != nfaces)
    {
        triangleEnergies->resize(nfaces);
        triangleEnergies->setZero();
    }

    for (int i = 0; i < nfaces; i++)
    {
        Matrix<scalar, 9, 1> dE;
        Matrix<scalar, 9, 9> hEexact;
        Matrix<scalar, 9, 9> hEinexact;
        Vector3i face = F.row(i);
        scalar senergy = triangleStretchingEnergy(V.row(face[0]), V.row(face[1]), V.row(face[2]),
            abars[i], 
            faceThicknesses[i], 
            params, 
            dEnergy ? &dE : NULL, 
            hEnergyExact ? &hEexact : NULL,
            hEnergyInexact ? &hEinexact : NULL
            );
        result += senergy;
        if (triangleEnergies)
            (*triangleEnergies)[i] += senergy;

        if (dEnergy)
        {
            for (int j = 0; j < 3; j++)
                dEnergy->segment<3>(3 * face[j]) += dE.segment<3>(3 * j);
        }
        if (hEnergyExact || hEnergyInexact)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (int x = 0; x < 3; x++)
                    {
                        for (int y = 0; y < 3; y++)
                        {
                            if(hEnergyExact)
                                hEnergyExact->push_back(Eigen::Triplet<scalar>(3 * face[j] + x, 3 * face[k] + y, hEexact(3*j+x, 3*k+y)));
                            if(hEnergyInexact)
                                hEnergyInexact->push_back(Eigen::Triplet<scalar>(3 * face[j] + x, 3 * face[k] + y, hEinexact(3 * j + x, 3 * k + y)));
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < nfaces; i++)
    {
        Matrix<scalar, 18, 1> dE;
        Matrix<scalar, 18, 18> hE;
        Vector3i face = F.row(i);
        Vector3i wing = faceWings.row(i);
        Matrix2m atemp, btemp;
   //     Matrix2m bbar;
   //     bbar.setZero();
        Vector3m wingverts[3];
        for (int j = 0; j < 3; j++)
            if (wing[j] != -1)
                wingverts[j] = V.row(wing[j]);
        scalar benergy = triangleBendingEnergy(V.row(face[0]), V.row(face[1]), V.row(face[2]),
            wing[0] == -1 ? NULL : &wingverts[0],
            wing[1] == -1 ? NULL : &wingverts[1],
            wing[2] == -1 ? NULL : &wingverts[2],
            abars[i], bbars[i], faceThicknesses[i], params, &atemp, &btemp, dEnergy ? &dE : NULL, (hEnergyExact || hEnergyInexact) ? &hE : NULL);        
        
        result += benergy;
        
        if (acurrent){
            acurrent->push_back(atemp);
        }
        if (bcurrent){
            bcurrent->push_back(btemp);
        }

        if (triangleEnergies)
            (*triangleEnergies)[i] += benergy;

        if (dEnergy)
        {
            for (int j = 0; j < 3; j++)
            {
                dEnergy->segment<3>(3 * face[j]) += dE.segment<3>(3 * j);
                if (wing[j] != -1)
                {
                    dEnergy->segment<3>(3 * wing[j]) += dE.segment<3>(9 + 3 * j);
                }
            }
        }
        if (hEnergyExact || hEnergyInexact)
        {
            for (int j = 0; j < 6; j++)
            {
                int jind;
                if (j < 3)
                    jind = 3 * face[j];
                else
                {
                    if (wing[j - 3] == -1)
                        continue;
                    jind = 3 * wing[j - 3];
                }
                for (int k = 0; k < 6; k++)
                {
                    int kind;
                    if (k < 3)
                        kind = 3 * face[k];
                    else
                    {
                        if (wing[k - 3] == -1)
                            continue;
                        kind = 3 * wing[k - 3];
                    }
                    for (int x = 0; x < 3; x++)
                    {
                        for (int y = 0; y < 3; y++)
                        {
                            if(hEnergyExact)
                                hEnergyExact->push_back(Eigen::Triplet<scalar>(jind + x, kind + y, hE(3 * j + x, 3 * k + y)));
                            if(hEnergyInexact)
                                hEnergyInexact->push_back(Eigen::Triplet<scalar>(jind + x, kind + y, hE(3 * j + x, 3 * k + y)));
                        }
                    }
                }
            }
        }
    }
    return result;
}

scalar poseEnergy(
    const MatrixX3m &V,
    const Eigen::MatrixX3i &F,
    const MatrixX3m &origV,
    VectorXm *dEnergy,
    std::vector<Eigen::Triplet<scalar> > *hEnergyExact,
    std::vector<Eigen::Triplet<scalar> > *hEnergyInexact)
{
    int nverts = (int)V.rows();

    if (dEnergy && dEnergy->size() != 3 * nverts)
    {
        dEnergy->resize(3 * nverts);
        dEnergy->setZero();
    }

    Vector3m comorig;
    comorig.setZero();
    Vector3m comcur;
    comcur.setZero();
    
    for (int i = 0; i < nverts; i++)
    {
        comorig += origV.row(i);
        comcur += V.row(i);
    }
    comorig /= nverts;
    comcur /= nverts;

    Vector3m ang;
    ang.setZero();

    Vector3m z(0, 0, 1);
    for (int i = 0; i < nverts; i++)
    {
        Vector3m origarm = origV.row(i).transpose() - comorig;
        Vector3m diff = V.row(i) - origV.row(i);
        ang += origarm.cross(diff);
    }
    ang /= nverts;

    scalar k = 1e8;
    scalar result = 0.5 * k * (comorig - comcur).squaredNorm();// +0.5*k*ang.squaredNorm();

    if (dEnergy)
    {
        for (int i = 0; i < nverts; i++)
        {
            dEnergy->segment<3>(3 * i) -= k * (comorig - comcur) / nverts;
            Vector3m origarm = origV.row(i).transpose() - comorig;
            //dEnergy->segment<3>(3 * i) -= k * origarm.cross(ang) / nverts;
        }
    }

    if (hEnergyExact || hEnergyInexact)
    {
        for (int i = 0; i < nverts; i++)
        {
            Vector3m origarm = origV.row(i).transpose() - comorig;
            Matrix3m cm = crossMatrix(origarm);
            Matrix3m Hang = -k*cm*cm / nverts / nverts;
            for (int j = 0; j < 3; j++)
            {
                if (hEnergyExact)
                    hEnergyExact->push_back(Eigen::Triplet<scalar>(3 * i + j, 3 * i + j, k / nverts / nverts));
                if (hEnergyInexact)
                    hEnergyInexact->push_back(Eigen::Triplet<scalar>(3 * i + j, 3 * i + j, k / nverts / nverts));
            }
            /*
            for (int j = 0; j < 3; j++)
            {
                for (int l = 0; l < 3; l++)
                {
                    if (hEnergyExact)
                        hEnergyExact->push_back(Eigen::Triplet<scalar>(3 * i + j, 3 * i + l, Hang(j,l)));
                    if (hEnergyInexact)
                        hEnergyInexact->push_back(Eigen::Triplet<scalar>(3 * i + j, 3 * i + l, Hang(j, l)));
                }
            }*/
        }
    }
    return result;
}
