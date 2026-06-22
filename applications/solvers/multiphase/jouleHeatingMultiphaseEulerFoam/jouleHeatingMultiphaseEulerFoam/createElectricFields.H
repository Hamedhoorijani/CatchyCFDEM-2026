/*---------------------------------------------------------------------------*\
    createElectricFields.H

    Field and property creation for the exact Kibret et al. electric model.
\*---------------------------------------------------------------------------*/

Info<< "Reading electricProperties\n" << endl;

IOdictionary electricProperties
(
    IOobject
    (
        "electricProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// -------------------------------------------------------------------------
// Model controls / names
// -------------------------------------------------------------------------
const word solidPhaseName
(
    electricProperties.lookupOrDefault<word>("solidPhaseName", "particles")
);

const word groundPatch
(
    electricProperties.lookupOrDefault<word>("groundPatch", "ground")
);

const word highPotentialPatch
(
    electricProperties.lookupOrDefault<word>("highPotentialPatch", "cathode")
);

const Switch useVoltageController
(
    electricProperties.lookupOrDefault<Switch>("useVoltageController", false)
);

// -------------------------------------------------------------------------
// Material / controller parameters
// -------------------------------------------------------------------------
const dimensionSet dimVoltage(1, 2, -3, 0, 0, -1, 0);
const dimensionSet dimConductivity(-1, -3, 3, 0, 0, 2, 0);
const dimensionSet dimEField(1, 1, -3, 0, 0, -1, 0);
const dimensionSet dimCurrentDensity(0, -2, 0, 0, 0, 1, 0);
const dimensionSet dimQJoule(1, -1, -3, 0, 0, 0, 0);
const dimensionSet dimVoltPerKelvin(1, 2, -3, -1, 0, -1, 0);

const dimensionedScalar sigmaParticle
(
    "sigmaParticle",
    dimConductivity,
    electricProperties.lookupOrDefault<scalar>("sigmaParticle", 25.0)
);

const dimensionedScalar desiredTemperature
(
    "desiredTemperature",
    dimTemperature,
    electricProperties.lookupOrDefault<scalar>("desiredTemperature", 1140.0)
);

const dimensionedScalar controllerGain
(
    "controllerGain",
    dimVoltPerKelvin,
    electricProperties.lookupOrDefault<scalar>("controllerGain", 0.5)
);

// -------------------------------------------------------------------------
// Electric potential field
// This should exist in 0/phiE with your boundary conditions.
// -------------------------------------------------------------------------
volScalarField phiE
(
    IOobject
    (
        "phiE",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// -------------------------------------------------------------------------
// Derived electric fields
// -------------------------------------------------------------------------
volScalarField sigmaEff
(
    IOobject
    (
        "sigmaEff",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroSigma", dimConductivity, 0.0)
);

volVectorField electricField
(
    IOobject
    (
        "electricField",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zeroE", dimEField, vector::zero)
);

volScalarField EMagnitude
(
    IOobject
    (
        "EMagnitude",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroEmag", dimEField, 0.0)
);

volVectorField currentDensity
(
    IOobject
    (
        "currentDensity",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zeroJ", dimCurrentDensity, vector::zero)
);

volScalarField JMagnitude
(
    IOobject
    (
        "JMagnitude",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroJmag", dimCurrentDensity, 0.0)
);

volScalarField QJoule
(
    IOobject
    (
        "QJoule",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroQJ", dimQJoule, 0.0)
);
