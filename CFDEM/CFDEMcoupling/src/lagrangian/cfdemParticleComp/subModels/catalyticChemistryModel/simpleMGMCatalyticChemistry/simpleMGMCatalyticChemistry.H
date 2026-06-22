#ifndef simpleMGMCatalyticChemistry_H
#define simpleMGMCatalyticChemistry_H

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "catalyticChemistryModel.H"
#include "cfdemCloudEnergy.H"
#include "rhoReactionThermo.H"
#include "basicGSChemistryModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "HashTable.H"
#include "List.H"
#include "word.H"

namespace Foam
{

class simpleMGMCatalyticChemistry : public catalyticChemistryModel
{
    //- Disallow default bitwise copy/assign
    simpleMGMCatalyticChemistry(const simpleMGMCatalyticChemistry&);
    void operator=(const simpleMGMCatalyticChemistry&);

private:

    // --- dictionaries / model handles
    const dictionary           propsDict_;
    const word                 catalystName_;

    autoPtr<rhoReactionThermo>     catThermoPtr_;
    autoPtr<basicGSChemistryModel> catChemistryPtr_;

    // --- fluid thermo refs
    const volScalarField& Tg_;
    const volScalarField& rhog_;
    const volScalarField& p_;

    // gas/solid composition fields
    const PtrList<volScalarField>& Yg_;   // gas mass-fractions from carrier thermo
    PtrList<volScalarField>        Ygs_;  // gas mass-fractions (catalyst thermo)
    PtrList<volScalarField>        Ys_;   // solid coverages (catalyst thermo)

    // species tables
    const speciesTable& solidSpecies_;
    const speciesTable& gasSpecies_;

    // sizes
    const label nGasSpecie_;
    const label nSolidSpecie_;

    // options
    const bool  useFluidTemperature_;
    const bool  useParticleQdot_;

    // particle/phase fields
    PtrList<volScalarField> RR_; // reaction rates for gas species
    volScalarField          Qdot_;

    // CFDEM exchange arrays (LIGGGHTS-facing)
    mutable double **partRR_;
    List<List<scalar> > ofPartRR_;
    
    mutable double **partCoverages_;
    List<List<scalar> > ofPartCoverages_;
    
    mutable double **partGas_;
    List<List<scalar> > ofPartGas_;
    
    mutable double **partPressure_;
    word     partTempName_;
    
    mutable double **partTemp_;
    
    mutable double **partQdot_;
    
    word     partHeatSourceName_;
    mutable double **partHeatSource_;
    
    word     partPressName_;

    bool     ini_;
    bool     massTransfer_;

    // particle diameter/mass IO
    word     partDName_;
    word     partMassName_;
    mutable double **partD_;
    mutable double **partMass_;
    scalar   growPelletExternally_;

    // transient buffers for chemistry
    scalarField dcdt_;
    scalarField Yi_;
    scalarField Ysi_;

    // --- MGM controls/state
    bool   mgmEnabled_;
    label  mgmEveryNSteps_;
    label  mgmSubCycles_;
    label  mgmNr_;
    label  nGrainLayers_;
    label  grainsPerLayer_;
    scalarField Rg0ByRp_;
    scalar deltaC0_;
    scalar porosity_;
    scalar tortuosity_;
    scalar kSolid_;
    scalar rhoSolid_;
    scalar CpSolid_;
    HashTable<scalar> DpCarbon_;
    HashTable<label>  nC_;
    scalar rhoCoke_;

    wordList  mgmSpeciesNames_;
    labelList mgmSpeciesIds_;

    List< List<scalar> > deltaC_;      // per-particle layer thicknesses
    scalarField          partRadius_;  // (optional) per-particle radius cache
    List< List<scalar> > mgmPoreN_;    // (optional) bookkeeping
    List< List<scalar> > mgmJsurf_;    // bookkeeping
    label                mgmStep_;
    List<scalar>         mgmBuf1_;
    List<scalar>         mgmBuf2_;

    // -- helpers
    void allocateMyArrays();

    void perLayerChemistryAndSinks_
    (
        const label pid,
        const scalar Ts,
        const List< List<scalar> >& CbyNode,
        const List<scalar>& aS_byCell,
        List< List<scalar> >& SbyNode,
        List<scalar>& qVolByCell
    );

    void solveEnergy1D_
    (
        const label pid,
        const scalar dt,
        const scalar TgCell,
        const scalar hgs,
        const List<scalar>& qVolByCell
    );

public:

    //- Runtime type information
    TypeName("simpleMGMCatalyticChemistry");

    //- Construct from components
    simpleMGMCatalyticChemistry(const dictionary& dict, cfdemCloudEnergy& sm);

    //- Destructor
    virtual ~simpleMGMCatalyticChemistry();

    //- Execute chemistry update
    virtual void execute();

    //- Post processing hook after flow (to push fields back to DEM)
    virtual void postFlow();

    //- Accessors
    const volScalarField RR(label i) const;
    const volScalarField Qdot() const;
};

} // End namespace Foam

#endif

