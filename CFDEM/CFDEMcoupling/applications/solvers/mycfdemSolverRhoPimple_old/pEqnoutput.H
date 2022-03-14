if ((!mesh.steady() && !pimple.simpleRho()) || pimple.consistent())
{
    rho = thermo.rho();
}

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

    Info<< "psip0 max/min: " << gMax(psip0) << "/" << gMin(psip0) << endl;

const volScalarField rAU("rAU", 1.0/UEqn.A());
    Info<< "rAU max/min: " << gMax(rAU) << "/" << gMin(rAU) << endl;
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(voidfraction*rho*rAU));

if (modelType=="A")
{
    rhorAUf *= fvc::interpolate(voidfraction);
}

    Info<< "rhorAUf max/min: " << gMax(rhorAUf) << "/" << gMin(rhorAUf) << endl;

tmp<volScalarField> rAtU
(
    pimple.consistent()
  ? volScalarField::New("rAtU", 1.0/(1.0/rAU - UEqn.H1()))
  : tmp<volScalarField>(nullptr)
);
tmp<surfaceScalarField> rhorAtUf
(
    pimple.consistent()
  ? surfaceScalarField::New("rhoRAtUf", fvc::interpolate(voidfraction*rho*rAtU()))
  : tmp<surfaceScalarField>(nullptr)
);

const volScalarField& rAAtU = pimple.consistent() ? rAtU() : rAU;
const surfaceScalarField& rhorAAtUf =
    pimple.consistent() ? rhorAtUf() : rhorAUf;
    Info<< "rhorAAtUf max/min: " << gMax(rhorAAtUf) << "/" << gMin(rhorAAtUf) << endl;

volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    Info<< "HbyA max/min: " << gMax(HbyA) << "/" << gMin(HbyA) << endl;

if (pimple.nCorrPiso() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiUs("phiUs", fvc::interpolate(voidfraction*rho*rAU*Ksl*Us)& mesh.Sf());
    Info<< "PhiUS max/min: " << gMax(phiUs) << "/" << gMin(phiUs) << endl;

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(voidfraction*rho)*fvc::flux(HbyA)
);
    Info<< "PhiHbyA max/min: " << gMax(phiHbyA) << "/" << gMin(phiHbyA) << endl;

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, U, phiHbyA, rhorAAtUf);

if (pimple.transonic())
{
    NotImplemented;
}
else
{
    if (pimple.consistent())
    {
        phiHbyA += (rhorAAtUf - rhorAUf)*fvc::snGrad(p)*mesh.magSf();
        HbyA += (rAAtU - rAU)*fvc::grad(p);
    }
    
    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(voidfraction, rho) + psi*correction(fvm::ddt(p))
      + fvc::div(phiHbyA)
     ==
        fvOptions(psi, p, rho.name())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAAtUf, p));

        pEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + phiUs + pEqn.flux();
        }
    }
}
    Info<< "Phi max/min: " << gMax(phi) << "/" << gMin(phi) << endl;

const bool limitedp = pressureControl.limit(p);

// Thermodynamic density update
thermo.correctRho(psi*p - psip0);

if (limitedp)
{
    rho = thermo.rho();
}

#include "rhoEqn.H"
    Info<< "Rho max/min: " << gMax(rho) << "/" << gMin(rho) << endl;
#include "compressibleContinuityErrsPU.H"

// Explicitly relax pressure for momentum corrector
p.relax();

if (modelType=="A")
{
    U = HbyA - rAU*(voidfraction*fvc::grad(p)-Ksl*Us);
}
else
{
    U = HbyA - rAU*(fvc::grad(p)-Ksl*Us);
}
    Info<< "U max/min: " << gMax(U) << "/" << gMin(U) << endl;
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(voidfraction, p);
}
    Info<< "P max/min: " << gMax(p) << "/" << gMin(p) << endl;
