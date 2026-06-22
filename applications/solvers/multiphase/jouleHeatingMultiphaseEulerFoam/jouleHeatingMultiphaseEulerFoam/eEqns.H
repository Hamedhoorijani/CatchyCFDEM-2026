const word solidPhaseName("particles");
const phaseModel& solid = fluid.phases()[solidPhaseName];

const volScalarField& alphaS = solid;
const volScalarField& Tsolid = solid.thermo().T();

tmp<volScalarField> trhoSolid = solid.rho();
const volScalarField& rhoSolid = trhoSolid();

// --------------------------------------------------------------------
// Exact paper conductivity:
//     sigmaEff = epsilon_s * sigma_p
// --------------------------------------------------------------------
sigmaEff = alphaS*sigmaParticle;
sigmaEff.correctBoundaryConditions();

// --------------------------------------------------------------------
// Optional: Eq. (22) voltage controller
//     V = n*(Tdesired - Tface)
// If you do not want the controller, remove this block and keep phiE BCs
// fixed in 0/phiE instead.
// --------------------------------------------------------------------
if (useVoltageController)
{
    const label cathodePatchID =
        mesh.boundaryMesh().findPatchID(highPotentialPatch);

    if (cathodePatchID >= 0)
    {
        scalarField TfaceField
        (
            Tsolid.boundaryField()[cathodePatchID].patchInternalField()
        );

        scalar TfaceSum = gSum(TfaceField);
        label nFaces = TfaceField.size();
        reduce(nFaces, sumOp<label>());

        scalar Tface =
            (nFaces > 0)
          ? TfaceSum/scalar(nFaces)
          : desiredTemperature.value();

        scalar Vnew =
            controllerGain.value()
           *(desiredTemperature.value() - Tface);

        refCast<fixedValueFvPatchScalarField>
        (
            phiE.boundaryFieldRef()[cathodePatchID]
        ) == Vnew;
    }
}

// Ground electrode at zero potential
{
    const label groundPatchID =
        mesh.boundaryMesh().findPatchID(groundPatch);

    if (groundPatchID >= 0)
    {
        refCast<fixedValueFvPatchScalarField>
        (
            phiE.boundaryFieldRef()[groundPatchID]
        ) == 0.0;
    }
}

phiE.correctBoundaryConditions();

// --------------------------------------------------------------------
// Solver-only regularization for OF8
// Without this, cells with alphaS = 0 give sigmaEff = 0 and the
// Laplacian matrix can have zero diagonal entries -> DIC/GAMG crashes.
// --------------------------------------------------------------------
const dimensionedScalar sigmaFloor
(
    "sigmaFloor",
    sigmaEff.dimensions(),
    electricProperties.lookupOrDefault<scalar>("sigmaFloor", 1e-12)
);

volScalarField sigmaSolve
(
    "sigmaSolve",
    max(sigmaEff, sigmaFloor)
);
sigmaSolve.correctBoundaryConditions();

Info<< "sigmaEff min/max = "
    << gMin(sigmaEff) << " " << gMax(sigmaEff) << nl
    << "sigmaSolve min/max = "
    << gMin(sigmaSolve) << " " << gMax(sigmaSolve) << endl;
// --------------------------------------------------------------------
// Exact paper potential equation:
//     div(sigmaEff * grad(phiE)) = 0
// --------------------------------------------------------------------
/*
fvScalarMatrix phiEEqn
(
    fvm::laplacian(sigmaEff, phiE)
);
phiEEqn.solve();
*/


solve(fvm::laplacian(sigmaSolve, phiE));

// --------------------------------------------------------------------
// Exact paper field/current/heating definitions:
//     E  = grad(phiE)
//     J  = sigmaEff * E
//     QJ = J . E
// --------------------------------------------------------------------
electricField = -fvc::grad(phiE);
EMagnitude    = mag(electricField);

currentDensity = sigmaEff*electricField;
JMagnitude     = mag(currentDensity);

QJoule = currentDensity & electricField;

// --------------------------------------------------------------------
// Diagnostics
// --------------------------------------------------------------------

surfaceScalarField sigmaEfff
(
    "sigmaEfff",
    fvc::interpolate(sigmaEff)
);

surfaceScalarField currentFlux
(
    "currentFlux",
    -sigmaEfff*fvc::snGrad(phiE)*mesh.magSf()
);

scalar I_highPot = 0.0;
scalar I_ground  = 0.0;
scalar V_applied = 0.0;

forAll(mesh.boundary(), patchI)
{
    const word& patchName = mesh.boundary()[patchI].name();
    const scalarField& patchFlux = currentFlux.boundaryField()[patchI];
    scalar patchCurrent = gSum(patchFlux);

    if (patchName == groundPatch)
    {
        I_ground += patchCurrent;
    }
    else if (patchName == highPotentialPatch)
    {
        I_highPot += patchCurrent;

        if (phiE.boundaryField()[patchI].size() > 0)
        {
            V_applied = gMax(phiE.boundaryField()[patchI]);
        }
    }
}


scalar I_total    = mag(I_ground);
scalar I_check    = mag(I_highPot);
scalar totalPower = fvc::domainIntegrate(QJoule).value();
scalar R_bed      = (I_total > SMALL) ? V_applied/I_total : GREAT;
scalar P_VI       = V_applied*I_total;

scalar totalSolidVol = fvc::domainIntegrate(alphaS).value();
scalar avgSigmaBed = 0.0;

if (totalSolidVol > SMALL)
{
    avgSigmaBed =
        fvc::domainIntegrate(sigmaEff).value()/totalSolidVol;
}

