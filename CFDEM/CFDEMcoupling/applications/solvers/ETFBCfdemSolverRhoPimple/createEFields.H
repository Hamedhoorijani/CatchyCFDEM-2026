Info<< "Reading electric potential field (V) \n" << endl;    
volScalarField V
(
    IOobject
    (
    "V",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
    ),
    mesh
);

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

volVectorField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimForce/dimCurrent/dimTime, Zero)
);

dimensionedScalar fluidElectricConductivity("fluidElectricConductivity", dimensionSet(-1, -3, 3, 0, 0, 2, 0),readScalar(electricProperties.lookup("fElecCond")));

