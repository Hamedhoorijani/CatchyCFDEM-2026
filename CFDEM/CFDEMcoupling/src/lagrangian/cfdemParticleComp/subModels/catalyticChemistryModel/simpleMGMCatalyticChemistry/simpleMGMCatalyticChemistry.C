/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "simpleMGMCatalyticChemistry.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "ListOps.H"      // Foam::findIndex
#include "zeroGradientFvPatchFields.H"

namespace Foam
{

defineTypeNameAndDebug(simpleMGMCatalyticChemistry, 0);
addToRunTimeSelectionTable(catalyticChemistryModel, simpleMGMCatalyticChemistry, dictionary);

// Small helper to clamp a label without invoking dimensioned overloads
static inline label clampLabel(label v, label lo, label hi)
{
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

// ----  MGM helpers --- //
static inline void assembleAreaDensityParticle
(
    const label mgmNr, const scalar Rp,
    const label nGrainLayers, const scalarField& Rg0ByRp,
    const List<scalar>& deltaCL, const label grainsPerLayer,
    List<scalar>& aS_byCell // size mgmNr
)
{
    const label N = (mgmNr > 1 ? mgmNr : 1);
    const scalar dr = Rp / scalar(N);
    forAll(aS_byCell, i) aS_byCell[i] = 0.0;

    for (label L=0; L<nGrainLayers; ++L)
    {
        const scalar Rg0 = Rg0ByRp[min(L,(label)Rg0ByRp.size()-1)] * Rp;
        const scalar Rc  = Rg0 + deltaCL[L];
        const scalar Ag  = 4.0*constant::mathematical::pi*sqr(Rg0)*grainsPerLayer; // per particle
        const scalar Vp  = (4.0/3.0)*constant::mathematical::pi*(Rp*Rp*Rp);
        const scalar aS  = Ag/max(SMALL,Vp); // [1/m]
        const label iCell = clampLabel(label(Rc/dr), 0, N-1);
        aS_byCell[iCell] += aS;
    }
}

// area density of one layer (used for Newton residual)
static inline scalar areaDensityOfLayer
(
    const scalar Rp, const scalar Rg0ByRpL, const label grainsPerLayer
)
{
    const scalar Rg0 = Rg0ByRpL * Rp;
    const scalar Ag  = 4.0*constant::mathematical::pi*sqr(Rg0)*grainsPerLayer;
    const scalar Vp  = (4.0/3.0)*constant::mathematical::pi*(Rp*Rp*Rp);
    return Ag/max(SMALL,Vp); // [1/m]
}

// diffusion resistance of product layer shell (spherical, exact)
static inline scalar shellResistance(const scalar Rg0, const scalar delta, const scalar Dp)
{
    if (delta <= 1e-12) return SMALL;
    const scalar R1 = Rg0, R2 = Rg0 + delta;
    return (R2 - R1)/max(SMALL, Dp*R1*R2); // [s/m]
}

// tiny dense solver (Newton): solve A x = b
static inline void solveDense(List< List<scalar> >& A, List<scalar>& b)
{
    const label n = b.size();
    for (label k=0; k<n; ++k)
    {
        scalar piv = A[k][k]; if (mag(piv) < SMALL) piv = SMALL;
        const scalar invp = 1.0/piv;
        for (label j=k; j<n; ++j) A[k][j] *= invp;
        b[k] *= invp;
        for (label i=0; i<n; ++i) if (i!=k)
        {
            const scalar f = A[i][k];
            for (label j=k; j<n; ++j) A[i][j] -= f*A[k][j];
            b[i] -= f*b[k];
        }
    }
}

static inline scalar pelletVolumeFromD(const scalar Dp)
{
    return (constant::mathematical::pi/6.0) * (Dp*Dp*Dp);
}

static inline scalar Deff_eps_tau(const scalar eps, const scalar tau, const scalar Dgas)
{
    return eps*max(SMALL, Dgas)/max(SMALL, tau);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

simpleMGMCatalyticChemistry::simpleMGMCatalyticChemistry
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    catalyticChemistryModel(dict, sm),
    propsDict_(dict.subDict(typeName + "Props")),
    catalystName_(propsDict_.lookup("catalystName")),
    catThermoPtr_(rhoReactionThermo::New(particleCloud_.mesh(), catalystName_)),
    catChemistryPtr_(basicGSChemistryModel::New(catThermoPtr_())),
    Tg_(sm.thermo().T()),
    rhog_(sm.thermo().rho()),
    p_(sm.thermo().p()),
    Yg_(dynamic_cast<const rhoReactionThermo&>(sm.thermo()).composition().Y()),
    Ygs_(catThermoPtr_->composition().Y()),
    Ys_(catThermoPtr_->composition().Ys()),
    solidSpecies_(catThermoPtr_->composition().solidSpecies()),
    gasSpecies_(catThermoPtr_->composition().species()),
    nGasSpecie_(Yg_.size()),
    nSolidSpecie_(Ys_.size()),
    useFluidTemperature_(propsDict_.lookupOrDefault("useFluidTemperature", false)),
    useParticleQdot_(propsDict_.lookupOrDefault("useParticleQdot", true)),
    RR_(nGasSpecie_),
    Qdot_
    (
        IOobject
        (
            IOobject::groupName("Qdot", catalystName_),
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        particleCloud_.mesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    partRR_(NULL),
    ofPartRR_(nGasSpecie_),
    partCoverages_(NULL),
    ofPartCoverages_(nSolidSpecie_),
    partGas_(NULL),
    ofPartGas_(nGasSpecie_),
    partPressure_(NULL),
    partTempName_(propsDict_.lookupOrDefault<word>("partTempName", "Temp")),
    partTemp_(NULL),
    partQdot_(NULL),
    partHeatSourceName_(propsDict_.lookupOrDefault<word>("partHeatSourceName", "heatSource")),
    partHeatSource_(NULL),
    partPressName_(propsDict_.lookupOrDefault<word>("partPressName", "p")),
    ini_(propsDict_.lookupOrDefault("initialized", false)),
    massTransfer_(false),
    partDName_(propsDict_.lookupOrDefault<word>("partDName", "d")),
    partMassName_(propsDict_.lookupOrDefault<word>("partMassName", "mass")),
    partD_(NULL),
    partMass_(NULL),
    growPelletExternally_(propsDict_.lookupOrDefault<scalar>("growPelletExternally", 0.0)),
    mgmEnabled_(false),
    mgmEveryNSteps_(20),
    mgmSubCycles_(1),
    mgmNr_(12),
    nGrainLayers_(3),
    grainsPerLayer_(8),
    Rg0ByRp_(nGrainLayers_, 0.2),
    deltaC0_(0.0),
    porosity_(0.35),
    tortuosity_(3.0),
    kSolid_(2.0),
    rhoSolid_(3000.0),
    CpSolid_(800.0),
    DpCarbon_(),
    nC_(),
    rhoCoke_(1800.0),
    mgmSpeciesNames_(),
    mgmSpeciesIds_(),
    deltaC_(),
    partRadius_(),
    mgmPoreN_(),
    mgmJsurf_(),
    mgmStep_(0),
    mgmBuf1_(),
    mgmBuf2_()
{
    if (sm.nrMassTransferModels() != 0) { massTransfer_ = true; }

    // allocate DEM arrays
    allocateMyArrays();

    // RR output fields
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("RR."+Yg_[fieldi].name(), catalystName_),
                    particleCloud_.mesh().time().timeName(),
                    particleCloud_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                particleCloud_.mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0),
                zeroGradientFvPatchField<scalar>::typeName
            )
        );
    }
    particleCloud_.checkCG(true);

    // MGM section
    if (propsDict_.found("MGM"))
    {
        const dictionary& mgmDict = propsDict_.subDict("MGM");

        mgmEnabled_     = mgmDict.lookupOrDefault<bool>("enableMGM", false);
        mgmEveryNSteps_ = mgmDict.lookupOrDefault<label>("mgmEveryNSteps", 20);
        mgmSubCycles_   = mgmDict.lookupOrDefault<label>("mgmSubCycles", 1);
        mgmNr_          = mgmDict.lookupOrDefault<label>("mgmNr", 12);

        if (mgmDict.found("speciesSubset")) mgmDict.lookup("speciesSubset") >> mgmSpeciesNames_;
        nGrainLayers_   = mgmDict.lookupOrDefault<label>("nGrainLayers", 3);
        grainsPerLayer_ = mgmDict.lookupOrDefault<label>("grainsPerLayer", 8);
        if (mgmDict.found("Rg0ByRp")) mgmDict.lookup("Rg0ByRp") >> Rg0ByRp_;
        deltaC0_        = mgmDict.lookupOrDefault<scalar>("deltaC0", 0.0);

        porosity_   = mgmDict.lookupOrDefault<scalar>("porosity", 0.35);
        tortuosity_ = mgmDict.lookupOrDefault<scalar>("tortuosity", 3.0);
        kSolid_     = mgmDict.lookupOrDefault<scalar>("kSolid", 2.0);
        rhoSolid_   = mgmDict.lookupOrDefault<scalar>("rhoSolid", 3000.0);
        CpSolid_    = mgmDict.lookupOrDefault<scalar>("CpSolid", 800.0);

        if (mgmDict.found("DpCarbon"))
        {
            const dictionary& Dc = mgmDict.subDict("DpCarbon");
            forAllConstIter(dictionary, Dc, iter)
            {
                DpCarbon_.insert(iter().keyword(), readScalar(iter().stream()));
            }
        }

        rhoCoke_ = mgmDict.lookupOrDefault<scalar>("rhoCoke", 1800.0);

        if (mgmDict.found("nC"))
        {
            const dictionary& nCdict = mgmDict.subDict("nC");
            forAllConstIter(dictionary, nCdict, iter)
            {
                nC_.insert(iter().keyword(), label(readScalar(iter().stream())));
            }
        }
    }

    // speciesSubset -> indices
    mgmSpeciesIds_.setSize(mgmSpeciesNames_.size());
    forAll(mgmSpeciesNames_, j)
    {
        //label id = Foam::findIndex(gasSpecies_, mgmSpeciesNames_[j]);
        label id = -1;
	for (label k = 0; k < gasSpecies_.size(); ++k)
	{
	    if (gasSpecies_[k] == mgmSpeciesNames_[j]) { id = k; break; }
	}
	if (id < 0)
	{
	    Warning<< "MGM: species '" << mgmSpeciesNames_[j]
		   << "' not present in gasSpecies_. It will be ignored.\n";
	}
	mgmSpeciesIds_[j] = id;

        
        if (id < 0)
        {
            Warning<< "MGM: species '" << mgmSpeciesNames_[j]
                   << "' not present in gasSpecies_. It will be ignored.\n";
        }
        mgmSpeciesIds_[j] = id;
    }

    // per-particle MGM state
    if (mgmEnabled_)
    {
        const label np = particleCloud_.numberOfParticles();
        deltaC_.setSize(np);
        partRadius_.setSize(np, 0.0);
        mgmPoreN_.setSize(np);
        mgmJsurf_.setSize(np);
        for (label p=0; p<np; ++p)
        {
            deltaC_[p].setSize(nGrainLayers_, deltaC0_);
            mgmPoreN_[p].setSize(mgmSpeciesIds_.size(), 0.0);
            mgmJsurf_[p].setSize(mgmSpeciesIds_.size(), 0.0);
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

simpleMGMCatalyticChemistry::~simpleMGMCatalyticChemistry()
{
    if (massTransfer_)
    {
        delete partGas_;
        delete partPressure_;
    }
    delete partTemp_;
    delete partQdot_;
    delete partRR_;
    delete partCoverages_;
    delete partMass_;
    delete partD_;
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * //

void simpleMGMCatalyticChemistry::allocateMyArrays()
{
    double initVal = 0.0;
    particleCloud_.dataExchangeM().allocateArray(partTemp_, initVal, 1);
    particleCloud_.dataExchangeM().allocateArray(partQdot_, initVal, 1);
    particleCloud_.dataExchangeM().allocateArray(partHeatSource_, initVal, 1);
    particleCloud_.dataExchangeM().allocateArray(partRR_, initVal, 1);
    particleCloud_.dataExchangeM().allocateArray(partCoverages_, initVal, 1);
    forAll(ofPartRR_, i) ofPartRR_[i].setSize(particleCloud_.numberOfParticles(), 0.0);
    forAll(ofPartCoverages_, i) ofPartCoverages_[i].setSize(particleCloud_.numberOfParticles(), Ys_[i][0]);

    if (massTransfer_)
    {
        particleCloud_.dataExchangeM().allocateArray(partPressure_, initVal, 1);
        forAll(ofPartGas_, i) ofPartGas_[i].setSize(particleCloud_.numberOfParticles(), Ygs_[i][0]);
        particleCloud_.dataExchangeM().allocateArray(partGas_, initVal, 1);
    };

    if (useParticleQdot_)
    {
        particleCloud_.dataExchangeM().getData(partHeatSourceName_, "scalar-atom", partHeatSource_);
    }

    particleCloud_.dataExchangeM().allocateArray(partD_,    initVal, 1);
    particleCloud_.dataExchangeM().allocateArray(partMass_, initVal, 1);

    // Read current diameter and mass from LIGGGHTS if available
    particleCloud_.dataExchangeM().getData(partDName_,    "scalar-atom", partD_);
    particleCloud_.dataExchangeM().getData(partMassName_, "scalar-atom", partMass_);
}


void simpleMGMCatalyticChemistry::perLayerChemistryAndSinks_
(
    const label pid, const scalar Ts,
    const List< List<scalar> >& CbyNode,
    const List<scalar>& aS_byCell,
    List< List<scalar> >& SbyNode,
    List<scalar>& qVolByCell
)
{
    // Reset outputs
    forAll(SbyNode, j) for (label i=0;i<mgmNr_;++i) SbyNode[j][i]=0.0;
    for (label i=0;i<mgmNr_;++i) qVolByCell[i]=0.0;

    const scalar Rp = 0.5*particleCloud_.d(pid);
    const label N = (mgmNr_ > 1 ? mgmNr_ : 1);
    const scalar dr = Rp / scalar(N);

    for (label L=0; L<nGrainLayers_; ++L)
    {
        // Locate macropore cell that hosts outer radius of this layer
        const scalar Rg0 = Rg0ByRp_[min(L,(label)Rg0ByRp_.size()-1)] * Rp;
        const scalar Rc  = Rg0 + deltaC_[pid][L];
        const label celli = clampLabel(label(Rc/dr), 0, N-1);
        const scalar aS   = areaDensityOfLayer(Rp, Rg0ByRp_[min(L,(label)Rg0ByRp_.size()-1)], grainsPerLayer_);

        // Build Cpore vector at this layer from macropore node
        List<scalar> Cpore(mgmSpeciesIds_.size(), 0.0);
        forAll(mgmSpeciesIds_, j) Cpore[j] = CbyNode[j][celli]; // [kmol/m3]

        // Unknowns: Cg (grain interface concentrations) for subset
        List<scalar> Cg = Cpore; // initial guess = Cpore

        // Newton iterations: enforce per-species balance Jdiff = Jrxn
        for (int it=0; it<8; ++it)
        {
            List<scalar> F(mgmSpeciesIds_.size(), 0.0);
            List< List<scalar> > J(mgmSpeciesIds_.size());
            forAll(J, r) J[r].setSize(mgmSpeciesIds_.size(), 0.0);

            forAll(mgmSpeciesIds_, j)
            {
                const label id = mgmSpeciesIds_[j]; if (id < 0) continue;
                const word sName = gasSpecies_[id];
                const scalar Dp = (DpCarbon_.found(sName) ? DpCarbon_[sName] : 1e-9); // [m2/s] diffusivity of species
                const scalar Rsh = shellResistance(Rg0, max(0.0, deltaC_[pid][L]), Dp); // [s/m]
                const scalar Jdiff = (Cpore[j] - Cg[j]) / max(SMALL, Rsh); // [kmol/m2/s]

                // chemistry kernel: Cg -> Yi_
                scalar Xden=0.0; 
                forAll(mgmSpeciesIds_, jj) Xden += Cg[jj];
                if (Xden <= SMALL) Xden = SMALL;
                forAll(Yi_, ii) Yi_[ii]=0.0; // clear
                forAll(mgmSpeciesIds_, jj)
                {
                    const label iid = mgmSpeciesIds_[jj]; if (iid < 0) continue;
                    const scalar Xjj = Cg[jj]/Xden;
                    const scalar Wij = catThermoPtr_->composition().Wi(iid);
                    Yi_[iid] = Xjj*Wij; // not normalized yet
                }
                scalar Ysum=VSMALL; 
                forAll(Yi_, ii) Ysum += Yi_[ii];
                forAll(Yi_, ii) Yi_[ii] /= Ysum;

                for (label ks=0; ks<nSolidSpecie_; ++ks) Ysi_[ks] = ofPartCoverages_[ks][pid];

                (void)catChemistryPtr_->getRatesQdotI
                (
                    particleCloud_.mesh().time().deltaTValue(),
                    particleCloud_.cellIDs()[pid][0],
                    dcdt_,
                    Yi_,
                    Ysi_,
                    rhog_[particleCloud_.cellIDs()[pid][0]],
                    Tg_[particleCloud_.cellIDs()[pid][0]],
                    Ts,
                    p_[particleCloud_.cellIDs()[pid][0]],
                    false
                );

                // volumetric sink -> surface flux
                List<scalar> Svol(mgmSpeciesIds_.size(), 0.0); // kmol/m3/s
                forAll(mgmSpeciesIds_, jj)
                {
                    const label iid = mgmSpeciesIds_[jj]; if (iid < 0) continue;
                    const scalar Wi = catThermoPtr_->composition().Wi(iid); // kg/kmol
                    Svol[jj] = - dcdt_[iid] / max(SMALL, Wi);
                }

                const scalar Jrxn = Svol[j] / max(SMALL, aS); // [kmol/m2/s]
                F[j] = Jdiff - Jrxn;
            }

            scalar nrm = 0.0; 
            forAll(F,r) nrm = max(nrm, mag(F[r]));
            if (nrm < 1e-6) break;

            // FD Jacobian
            forAll(Cg, jcol)
            {
                const scalar base = Cg[jcol];
                const scalar h = max(1e-8, 1e-3*max(1.0, base));
                Cg[jcol] = max(0.0, base + h);

                List<scalar> Fp(mgmSpeciesIds_.size(), 0.0);
                forAll(mgmSpeciesIds_, j)
                {
                    const label id = mgmSpeciesIds_[j]; if (id < 0) continue;
                    const word sName = gasSpecies_[id];
                    const scalar Dp = (DpCarbon_.found(sName) ? DpCarbon_[sName] : 1e-9);
                    const scalar Rsh = shellResistance(Rg0, max(0.0, deltaC_[pid][L]), Dp);
                    const scalar Jdiff = (Cpore[j] - Cg[j]) / max(SMALL, Rsh);

                    scalar Xden=0.0; 
                    forAll(mgmSpeciesIds_, jj) Xden += Cg[jj];
                    if (Xden <= SMALL) Xden = SMALL;
                    forAll(Yi_, ii) Yi_[ii]=0.0;
                    forAll(mgmSpeciesIds_, jj)
                    {
                        const label iid = mgmSpeciesIds_[jj]; if (iid < 0) continue;
                        const scalar Xjj = Cg[jj]/Xden;
                        const scalar Wij = catThermoPtr_->composition().Wi(iid);
                        Yi_[iid] = Xjj*Wij;
                    }
                    scalar Ysum=VSMALL; forAll(Yi_, ii) Ysum += Yi_[ii];
                    forAll(Yi_, ii) Yi_[ii] /= Ysum;
                    for (label ks=0; ks<nSolidSpecie_; ++ks) Ysi_[ks] = ofPartCoverages_[ks][pid];

                    (void)catChemistryPtr_->getRatesQdotI
                    (
                        particleCloud_.mesh().time().deltaTValue(),
                        particleCloud_.cellIDs()[pid][0],
                        dcdt_,
                        Yi_,
                        Ysi_,
                        rhog_[particleCloud_.cellIDs()[pid][0]],
                        Tg_[particleCloud_.cellIDs()[pid][0]],
                        Ts,
                        p_[particleCloud_.cellIDs()[pid][0]],
                        false
                    );

                    List<scalar> Svol(mgmSpeciesIds_.size(), 0.0);
                    forAll(mgmSpeciesIds_, jj)
                    {
                        const label iid = mgmSpeciesIds_[jj]; if (iid < 0) continue;
                        const scalar Wi = catThermoPtr_->composition().Wi(iid);
                        Svol[jj] = - dcdt_[iid] / max(SMALL, Wi);
                    }

                    const scalar Jrxn = Svol[j] / max(SMALL, aS);
                    Fp[j] = Jdiff - Jrxn;
                }
                Cg[jcol] = base;

                forAll(F, r) J[r][jcol] = (Fp[r]-F[r])/h;
            }

            List<scalar> dx = -F;
            solveDense(J, dx);
            forAll(Cg, j) Cg[j] = max(0.0, Cg[j] + dx[j]);
        } // Newton

        forAll(mgmSpeciesIds_, j)
        {
            const label id = mgmSpeciesIds_[j]; if (id < 0) continue;
            const word sName = gasSpecies_[id];
            const scalar Dp = (DpCarbon_.found(sName) ? DpCarbon_[sName] : 1e-9);
            const scalar Rsh = shellResistance(Rg0, max(0.0, deltaC_[pid][L]), Dp);
            const scalar Jdiff = (CbyNode[j][celli] - Cg[j]) / max(SMALL, Rsh); // [kmol/m2/s]
            SbyNode[j][celli] += aS * Jdiff; // [kmol/m3/s]
        }

        scalar Xden=0.0; 
        forAll(mgmSpeciesIds_, jj) Xden += Cg[jj];
        if (Xden <= SMALL) Xden = SMALL;
        forAll(Yi_, ii) Yi_[ii]=0.0;
        forAll(mgmSpeciesIds_, jj)
        {
            const label iid = mgmSpeciesIds_[jj]; if (iid < 0) continue;
            const scalar Xjj = Cg[jj]/Xden;
            const scalar Wij = catThermoPtr_->composition().Wi(iid);
            Yi_[iid] = Xjj*Wij;
        }
        scalar Ysum=VSMALL; 
        forAll(Yi_, ii) Ysum += Yi_[ii];
        forAll(Yi_, ii) Yi_[ii] /= Ysum;
        for (label ks=0; ks<nSolidSpecie_; ++ks) Ysi_[ks] = ofPartCoverages_[ks][pid];

        const scalar Qvol =
            catChemistryPtr_->getRatesQdotI
            (
                particleCloud_.mesh().time().deltaTValue(),
                particleCloud_.cellIDs()[pid][0],
                dcdt_,
                Yi_,
                Ysi_,
                rhog_[particleCloud_.cellIDs()[pid][0]],
                Tg_[particleCloud_.cellIDs()[pid][0]],
                Ts,
                p_[particleCloud_.cellIDs()[pid][0]],
                false
            );

        qVolByCell[celli] += Qvol; // [W/m3]
    } // layers
}


//Simplified radial energy equation - 1D

void simpleMGMCatalyticChemistry::solveEnergy1D_
(
    const label pid, 
    const scalar dt,
    const scalar TgCell, 
    const scalar hgs,
    const List<scalar>& qVolByCell
)
{
    const scalar Rp = 0.5*particleCloud_.d(pid);	//initial particle radius (pre-MGM @ the current step)
    const label N = (mgmNr_ > 1 ? mgmNr_ : 1);		// number of layers
    const scalar dr = Rp / scalar(N);			//discretization
    
    // initial profile = last avg or Tg
    List<scalar> T(N, partTemp_[pid][0] > 0 ? partTemp_[pid][0] : TgCell);

    // Build tridiagonal for implicit conduction with volumetric heat
    List<scalar> a(N,0.0), b(N,0.0), c(N,0.0), rhs(N,0.0);
    for (label i=0;i<N;++i)
    {
        const scalar rL=i*dr, rR=(i+1)*dr;
        const scalar AL=4.0*constant::mathematical::pi*sqr(rL);
        const scalar AR=4.0*constant::mathematical::pi*sqr(rR);
        const scalar VL=(4.0/3.0)*constant::mathematical::pi*((rR*rR*rR)-(rL*rL*rL));

        scalar aW=0.0, aE=0.0;
        if (i>0)   aW = kSolid_*AL/dr;
        if (i<N-1) aE = kSolid_*AR/dr;

        a[i]   = -aW;
        b[i]   = rhoSolid_*CpSolid_*VL/dt + aW + aE;
        c[i]   = -aE;
        rhs[i] = rhoSolid_*CpSolid_*VL/dt*T[i] + qVolByCell[i]*VL;

        if (i==0) a[i]=0.0; // symmetry
    }

    // Robin at surface: -k dT/dr = hgs (T - Tg)
    {
        const label i = N-1;
        const scalar rR=(i+1)*dr;
        const scalar AR=4.0*constant::mathematical::pi*sqr(rR);
        b[i]   += hgs*AR;
        rhs[i] += hgs*AR*TgCell;
    }

    // Solve tri-di
    for (label i=1;i<N;++i)
    {
        const scalar m = a[i]/max(SMALL,b[i-1]);
        b[i]   -= m*c[i-1];
        rhs[i] -= m*rhs[i-1];
    }
    T[N-1] = rhs[N-1]/max(SMALL,b[N-1]);
    for (label i=N-2;i>=0;--i) T[i] = (rhs[i]-c[i]*T[i+1])/max(SMALL,b[i]);

    // Cache volume-avg temperature into partTemp_ (LIGGGHTS-facing array)
    scalar Tav=0.0, Vsum=0.0;
    for (label i=0;i<N;++i)
    {
        const scalar rL=i*dr, rR=(i+1)*dr;
        const scalar VL=(4.0/3.0)*constant::mathematical::pi*((rR*rR*rR)-(rL*rL*rL));
        Tav  += T[i]*VL; Vsum += VL;
    }
    partTemp_[pid][0] = Tav/max(SMALL,Vsum);
}


void simpleMGMCatalyticChemistry::execute()
{
    forAll(Ygs_, i)
    {
        forAll(Ygs_[i], celli)
        {
            Ygs_[i][celli] = Yg_[i][celli];
        }
    }

    // initialize surface
    if (!ini_)
    {
        allocateMyArrays();

        label celli = 0;
        scalar pi = p_[celli];
        scalar Tgi = Tg_[celli];
        scalar Tsi = catThermoPtr_->T()[celli];

        forAll(Yi_, i) Yi_[i] = Ygs_[i][celli];
        forAll(Ysi_, i) Ysi_[i] = Ys_[i][celli];

        Qdot_[celli] =
            catChemistryPtr_->getRatesQdotI
            (
                particleCloud_.mesh().time().deltaTValue(),
                celli,
                dcdt_,
                Yi_,
                Ysi_,
                rhog_[celli],
                Tgi,
                Tsi,
                pi,
                true,  // initialize coverages
                !useParticleQdot_
            );

        for (int index = 0; index < particleCloud_.numberOfParticles(); ++index)
        {
            label cellI = particleCloud_.cellIDs()[index][0];
            if (cellI >= 0)
            {
                for (label i=0; i<nSolidSpecie_; i++) ofPartCoverages_[i][index] = Ysi_[i];
            }
        }

        if (massTransfer_)
        {
            forAll(ofPartGas_, i)
            {
                particleCloud_.dataExchangeM().getData(gasSpecies_[i], "scalar-atom", partGas_);
                forAll(ofPartGas_[i], index)
                {
                    ofPartGas_[i][index] = partGas_[index][0];
                }
            }
            particleCloud_.dataExchangeM().getData(partPressName_, "scalar-atom", partPressure_);
        }

        ini_ = true;
    }
    else
    {
        allocateMyArrays();

        forAll(ofPartCoverages_, i)
        {
            particleCloud_.dataExchangeM().getData(solidSpecies_[i], "scalar-atom", partCoverages_);
            forAll(ofPartCoverages_[i], index)
            {
                ofPartCoverages_[i][index] = partCoverages_[index][0];
            }
        }

        if (massTransfer_)
        {
            forAll(ofPartGas_, i)
            {
                particleCloud_.dataExchangeM().getData(gasSpecies_[i], "scalar-atom", partGas_);
                forAll(ofPartGas_[i], index)
                {
                    ofPartGas_[i][index] = partGas_[index][0];
                }
            }
            particleCloud_.dataExchangeM().getData(partPressName_, "scalar-atom", partPressure_);
        };

        if (!useFluidTemperature_)
        {
            particleCloud_.dataExchangeM().getData(partTempName_, "scalar-atom", partTemp_);
        }
    }

    for (int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        label cellI = particleCloud_.cellIDs()[index][0];
        if (cellI >= 0)
        {
            scalar pi;
            if (!massTransfer_)
            {
                for (label i=0; i<nGasSpecie_; i++) Yi_[i] = Ygs_[i][cellI];
                pi = p_[cellI];
            }
            else
            {
                for (label i=0; i<nGasSpecie_; i++) Yi_[i] = ofPartGas_[i][index];
                pi = partPressure_[index][0];
            }
            for (label i=0; i<nSolidSpecie_; i++) Ysi_[i] = ofPartCoverages_[i][index];

            if (useFluidTemperature_)
            {
                scalar Tgi = Tg_[cellI];

                partQdot_[index][0] =
                    catChemistryPtr_->getRatesQdotI
                    (
                        particleCloud_.mesh().time().deltaTValue(),
                        cellI,
                        dcdt_,
                        Yi_,
                        Ysi_,
                        rhog_[cellI],
                        Tg_[cellI],
                        Tgi,
                        pi,
                        false
                    );
                partQdot_[index][0] *= particleCloud_.particleVolume(index);
            }
            else
            {
                partQdot_[index][0] =
                    catChemistryPtr_->getRatesQdotI
                    (
                        particleCloud_.mesh().time().deltaTValue(),
                        cellI,
                        dcdt_,
                        Yi_,
                        Ysi_,
                        rhog_[cellI],
                        Tg_[cellI],
                        partTemp_[index][0],
                        pi,
                        false,
                        !useParticleQdot_
                    );
                partQdot_[index][0] *= particleCloud_.particleVolume(index);
            }
            for (label i=0; i<nGasSpecie_; i++) ofPartRR_[i][index] = dcdt_[i]*particleCloud_.particleVolume(index);
            for (label i=0; i<nSolidSpecie_; i++) ofPartCoverages_[i][index] = Ysi_[i];
        }
    }

    forAll(Ys_, i)
    {
        forAll(ofPartCoverages_[i], index)
        {
            partCoverages_[index][0] = ofPartCoverages_[i][index];
        }
        Ys_[i].primitiveFieldRef() = 0.0;
        particleCloud_.averagingM().resetWeightFields();
        particleCloud_.averagingM().setScalarAverage
        (
            Ys_[i],
            partCoverages_,
            particleCloud_.particleWeights(),
            particleCloud_.averagingM().UsWeightField(),
            NULL
        );
        Info << Ys_[i].name() << " max/min: " << max(Ys_[i]).value() << "/" << min(Ys_[i]).value() << endl;
    }

    Qdot_.primitiveFieldRef() = 0.0;
    if (!useParticleQdot_)
    {
        particleCloud_.averagingM().setScalarSum
        (
            Qdot_,
            partQdot_,
            particleCloud_.particleWeights(),
            NULL
        );
        Qdot_.primitiveFieldRef() /= Qdot_.mesh().V();
    }

    if (!massTransfer_)
    {
        forAll(RR_, i)
        {
            forAll(ofPartRR_[i], index)
            {
                partRR_[index][0] = ofPartRR_[i][index];
            }
            RR_[i].primitiveFieldRef() = 0.0;
            particleCloud_.averagingM().setScalarSum
            (
                RR_[i],
                partRR_,
                particleCloud_.particleWeights(),
                NULL
            );
            RR_[i].primitiveFieldRef() /= RR_[i].mesh().V();
        }
    }
    else
    {
        forAll(RR_, i) { RR_[i].primitiveFieldRef() = 0.0; }

        for (int index = 0; index < particleCloud_.numberOfParticles(); ++index)
        {
            label cellI = particleCloud_.cellIDs()[index][0];
            if (cellI >= 0)
            {
                const scalar dt = particleCloud_.mesh().time().deltaTValue();

                if (mgmEnabled_ && mgmSpeciesIds_.size())
                {
                    // Build bulk concentrations [kmol/m3] for subset from gas cell
                    List<scalar> Cbulk(mgmSpeciesIds_.size(), 0.0);
                    {
                        scalar YmW = 0.0;
                        forAll(Yg_, i) YmW += Yg_[i][cellI]/catThermoPtr_->composition().Wi(i);
                        forAll(mgmSpeciesIds_, j)
                        {
                            const label id = mgmSpeciesIds_[j];
                            if (id < 0) continue;
                            const scalar Yi = Yg_[id][cellI];
                            const scalar Wi = catThermoPtr_->composition().Wi(id);
                            const scalar Xi = Yi/Wi / max(SMALL, YmW);
                            const scalar Ct = rhog_[cellI] / max(SMALL, (YmW>0? 1.0/YmW : GREAT));
                            Cbulk[j] = Ct*Xi;
                        }
                    }

                    // Film coefficients for subset: fallback to large (Dirichlet at surface)
                    List<scalar> kgSubset(mgmSpeciesIds_.size(), GREAT);

                    const label N  = (mgmNr_ > 1 ? mgmNr_ : 1);
                    const scalar Rp = 0.5*particleCloud_.d(index);
                    const scalar dr = Rp/scalar(N);

                    // Estimate h_gs for this particle (uses last-step convective magnitude)
                    scalar hgs = 1e3; // [W/m2/K] fallback
                    {
                        const scalar d = particleCloud_.d(index);
                        const scalar As = constant::mathematical::pi * (d*d);
                        const scalar Tp = useFluidTemperature_ ? Tg_[cellI] : partTemp_[index][0];
                        const scalar dT = max(1e-3, mag(Tg_[cellI] - Tp));
                        const scalar qconv = max(0.0, mag(partHeatSource_[index][0]));
                        if (As > SMALL) hgs = qconv/(As*dT);
                    }

                    // Macropore fields
                    List< List<scalar> > C(mgmSpeciesIds_.size(), List<scalar>(N, 0.0));
                    List< List<scalar> > S(mgmSpeciesIds_.size(), List<scalar>(N, 0.0));
                    forAll(mgmSpeciesIds_, j) for (label iNode=0;iNode<N;++iNode) C[j][iNode]=Cbulk[j];

                    // Internal surface area density per macropore cell [1/m]
                    List<scalar> aS_byCell(N, 0.0);
                    assembleAreaDensityParticle(N, Rp, nGrainLayers_, Rg0ByRp_, deltaC_[index], grainsPerLayer_, aS_byCell);

                    // Heat source density per macropore cell [W/m3]
                    List<scalar> qVolByCell(N, 0.0);

                    // Picard coupling iterations
                    for (int it=0; it<3; ++it)
                    {
                        const scalar Tloc = useFluidTemperature_ ? Tg_[cellI] : partTemp_[index][0];
                        perLayerChemistryAndSinks_(index, Tloc, C, aS_byCell, S, qVolByCell);

                        // --- Coke growth from carbon-atom sinks
                        scalar dV_coke_total = 0.0;
                        for (label L=0; L<nGrainLayers_; ++L)
                        {
                            const scalar Rg0   = Rg0ByRp_[min(L,(label)Rg0ByRp_.size()-1)] * Rp;
                            const scalar Rc    = Rg0 + deltaC_[index][L];
                            const label celliL = clampLabel(label(Rc/dr), 0, N-1);

                            const scalar aS_L = areaDensityOfLayer(Rp, Rg0ByRp_[min(L,(label)Rg0ByRp_.size()-1)], grainsPerLayer_);

                            scalar S_C_vol = 0.0;
                            forAll(mgmSpeciesIds_, j)
                            {
                                const label gid = mgmSpeciesIds_[j];
                                if (gid < 0) continue;
                                const word sName = gasSpecies_[gid];
                                const label nC = nC_.found(sName) ? nC_[sName] : 0;
                                S_C_vol += nC * S[j][celliL];
                            }

                            if (S_C_vol > SMALL && aS_L > SMALL)
                            {
                                const scalar J_C = S_C_vol / aS_L;      // [kmol/m2/s]
                                const scalar dDelta = (J_C * 12.0 / rhoCoke_) * dt; // [m]

                                const scalar Rg0L = Rg0ByRp_[min(L,(label)Rg0ByRp_.size()-1)] * Rp;
                                const scalar A_layer = 4.0*constant::mathematical::pi*sqr(Rg0L) * grainsPerLayer_;
                                const scalar dV_layer = A_layer * dDelta;

                                dV_coke_total += max(0.0, dV_layer);

                                const scalar maxDelta = max(0.0, Rp - Rg0L - 1e-9);
                                deltaC_[index][L] = min(deltaC_[index][L] + dDelta, maxDelta);
                            }
                        }

                        if (it == 2 && dV_coke_total > 0)
                        {
                            partMass_[index][0] += rhoCoke_ * dV_coke_total;

                            if (growPelletExternally_ > SMALL)
                            {
                                const scalar dOld = partD_[index][0]; //pre-MGM diameter
                                const scalar Vold = pelletVolumeFromD(dOld); // pre-MGM particle volume
                                const scalar Vnew = Vold + growPelletExternally_ * dV_coke_total; 	//post-MGM particle diameter
                                // invert π/6 d^3
                                partD_[index][0]  = cbrt(6.0*Vnew/constant::mathematical::pi); 		//update the new particle diameter
                            }
                        }

                        forAll(aS_byCell, iA) aS_byCell[iA] = 0.0;
                        assembleAreaDensityParticle(N, Rp, nGrainLayers_, Rg0ByRp_, deltaC_[index], grainsPerLayer_, aS_byCell);

                        // 2) macropore implicit diffusion with sinks S[j]
                        forAll(mgmSpeciesIds_, j)
                        {
                            const scalar Dgas = 1e-5; // [m2/s] placeholder
                            const scalar De   = Deff_eps_tau(porosity_, tortuosity_, Dgas); // [m2/s]
                            const scalar kg   = kgSubset[j]; // [m/s]

                            List<scalar> a(N,0.0), b(N,0.0), c_(N,0.0), rhs(N,0.0);
                            for (label iNode=0;iNode<N;++iNode)
                            {
                                const scalar rL=iNode*dr, rR=(iNode+1)*dr;
                                const scalar AL=4.0*constant::mathematical::pi*sqr(rL);
                                const scalar AR=4.0*constant::mathematical::pi*sqr(rR);
                                const scalar VL=(4.0/3.0)*constant::mathematical::pi*((rR*rR*rR)-(rL*rL*rL));

                                scalar aW=0.0, aE=0.0;
                                if (iNode>0)   aW = De*AL/dr;
                                if (iNode<N-1) aE = De*AR/dr;

                                a[iNode]   = -aW;
                                b[iNode]   = porosity_*VL/dt + aW + aE;
                                c_[iNode]  = -aE;
                                rhs[iNode] = porosity_*VL/dt*C[j][iNode] - S[j][iNode]*VL;

                                if (iNode==0) a[iNode]=0.0; // symmetry
                            }
                            // Robin at surface
                            {
                                const label i = N-1;
                                const scalar rR=(i+1)*dr;
                                const scalar AR=4.0*constant::mathematical::pi*sqr(rR);
                                b[i]   += kg*AR;
                                rhs[i] += kg*AR*Cbulk[j];
                            }
                            // solve tri-di
                            for (label i=1;i<N;++i)
                            {
                                const scalar m = a[i]/max(SMALL,b[i-1]);
                                b[i]   -= m*c_[i-1];
                                rhs[i] -= m*rhs[i-1];
                            }
                            C[j][N-1] = rhs[N-1]/max(SMALL,b[N-1]);
                            for (label i=N-2;i>=0;--i) C[j][i] = (rhs[i]-c_[i]*C[j][i+1])/max(SMALL,b[i]);
                        }
                    } // Picard

                    // Update particle pore composition mass-fractions from MGM average - Simplified for unresolved CFD-DEM
                    List<scalar> nSubset(mgmSpeciesIds_.size(), 0.0);
                    scalar VporeTot = 0.0;
                    for (label iNode=0;iNode<N;++iNode)
                    {
                        const scalar rL=iNode*dr, rR=(iNode+1)*dr;
                        const scalar VL=(4.0/3.0)*constant::mathematical::pi*((rR*rR*rR)-(rL*rL*rL));
                        const scalar VporeCell = porosity_*VL;
                        VporeTot += VporeCell;
                        forAll(mgmSpeciesIds_, j) nSubset[j] += C[j][iNode]*VporeCell;
                    }

                    List<scalar> Xfull(nGasSpecie_, 0.0);
                    scalar denom = 0.0;
                    forAll(Yg_, i)
                    {
                        if (Foam::findIndex(mgmSpeciesIds_, i) >= 0) continue;
                        const scalar Wi = catThermoPtr_->composition().Wi(i);
                        Xfull[i] = Yg_[i][cellI]/Wi;
                        denom += Xfull[i];
                    }
                    scalar nTot = 0.0; 
                    forAll(nSubset, j) nTot += nSubset[j];
                    forAll(mgmSpeciesIds_, j)
                    {
                        const label id = mgmSpeciesIds_[j];
                        if (id < 0) continue;
                        Xfull[id] = nSubset[j]/max(SMALL, nTot);
                    }
                    denom = 0.0; 
                    forAll(Xfull, i) denom += Xfull[i];
                    forAll(Xfull, i) Xfull[i] /= max(SMALL, denom);

                    scalar Mt = 0.0; 
                    forAll(Xfull, i) Mt += Xfull[i]*catThermoPtr_->composition().Wi(i);
                    forAll(Yg_, i) ofPartGas_[i][index] = Xfull[i]*catThermoPtr_->composition().Wi(i)/max(SMALL, Mt);

                    // Update particle pressure (simple ideal-gas balance)
                    scalar nDotTot = 0.0;
                    forAll(mgmSpeciesIds_, j)
                    {
                        for (label iNode=0;iNode<N;++iNode)
                        {
                            const scalar rL=iNode*dr, rR=(iNode+1)*dr;
                            const scalar VL=(4.0/3.0)*constant::mathematical::pi*((rR*rR*rR)-(rL*rL*rL));
                            const scalar VporeCell = porosity_*VL;
                            nDotTot += - S[j][iNode] * VporeCell;
                        }
                    }
                    const scalar Rbar = 8.314462618e3; // [kPa m3/(kmol K)]
                    scalar Tloc = useFluidTemperature_ ? Tg_[cellI] : partTemp_[index][0];
                    if (Tloc <= SMALL) Tloc = Tg_[cellI];
                    partPressure_[index][0] += (nDotTot * Rbar * Tloc * dt)/ max(SMALL, VporeTot);

                    // Intraparticle conduction with volumetric heat source qVolByCell
                    solveEnergy1D_(index, dt, Tg_[cellI], hgs, qVolByCell);
                }
                else
                {
                    // Original lumped update (fallback)
                    scalar RRnet_ = 0.0;
                    scalar R_ = 8.314462618;
                    scalar T_;
                    scalar Y_MW = 0;
                    if (useFluidTemperature_) { T_ = Tg_[cellI]; }
                    else { T_ = partTemp_[index][0]; }

                    forAll(Yg_, i) { Y_MW += ofPartGas_[i][index]/catThermoPtr_->composition().Wi(i); }
                    forAll(Yg_, i) { ofPartGas_[i][index] = ofPartGas_[i][index]/catThermoPtr_->composition().Wi(i)/Y_MW; }

                    forAll(Yg_, i)
                    {
                        RRnet_  += ofPartRR_[i][index]*dt/(catThermoPtr_->composition().Wi(i)*1e-3);
                        ofPartGas_[i][index] = max(ofPartGas_[i][index]*partPressure_[index][0]/(R_*T_)*catChemistryPtr_->porosity()*particleCloud_.particleVolume(index)
                                                   + ofPartRR_[i][index]*dt/(catThermoPtr_->composition().Wi(i)*1e-3), 0.0);
                    }

                    scalar Mt_ = 0.0;
                    scalar nt_ = 0.0;
                    forAll(Yg_, i) { Mt_ += ofPartGas_[i][index]*catThermoPtr_->composition().Wi(i); nt_ += ofPartGas_[i][index]; }
                    forAll(Yg_, i) { ofPartGas_[i][index] = ofPartGas_[i][index]*catThermoPtr_->composition().Wi(i)/Mt_; }

                    partPressure_[index][0] += RRnet_*8.314462618*T_/(catChemistryPtr_->porosity()*particleCloud_.particleVolume(index));
                    if (index==0) Info << "partPressure_: " << partPressure_[index][0] << endl;
                }
            }
            else
            {
                partPressure_[index][0] = 0;
            }
        }

        forAll(ofPartGas_, i)
        {
            forAll(ofPartGas_[i], index)
            {
                partGas_[index][0] = particleCloud_.cellIDs()[index][0] >= 0 ? ofPartGas_[i][index] : 0.0;
            }
            particleCloud_.dataExchangeM().giveData(gasSpecies_[i], "scalar-atom", partGas_);
        }
        particleCloud_.dataExchangeM().giveData(partPressName_, "scalar-atom", partPressure_);
        particleCloud_.dataExchangeM().giveData(partMassName_, "scalar-atom", partMass_);
        particleCloud_.dataExchangeM().giveData(partDName_,    "scalar-atom", partD_); 	//comunicate back the data to LIGGGHTS ( question is if LIGGGHTS can understand the data or not?)
    }
}

const volScalarField simpleMGMCatalyticChemistry::RR(label i) const
{
    return RR_[i];
}

const volScalarField simpleMGMCatalyticChemistry::Qdot() const
{
    return Qdot_;
}

void simpleMGMCatalyticChemistry::postFlow()
{
    forAll(ofPartCoverages_, i)
    {
        forAll(ofPartCoverages_[i], index)
        {
            partCoverages_[index][0] = particleCloud_.cellIDs()[index][0] >= 0 ? ofPartCoverages_[i][index] : 0.0;
        }
        particleCloud_.dataExchangeM().giveData(solidSpecies_[i], "scalar-atom", partCoverages_);
    }

    if (useParticleQdot_)
    {
        particleCloud_.dataExchangeM().getData(partHeatSourceName_,"scalar-atom", partHeatSource_);
        for (int index = 0; index < particleCloud_.numberOfParticles(); ++index)
        {
            label cellI = particleCloud_.cellIDs()[index][0];
            if (cellI >= 0)
            {
                partHeatSource_[index][0] = partHeatSource_[index][0] + partQdot_[index][0];
            }
            else
            {
                partHeatSource_[index][0] = 0.0;
            }
        }
        particleCloud_.dataExchangeM().giveData(partHeatSourceName_, "scalar-atom", partHeatSource_);
    }
}

// ************************************************************************* //

} // End namespace Foam

