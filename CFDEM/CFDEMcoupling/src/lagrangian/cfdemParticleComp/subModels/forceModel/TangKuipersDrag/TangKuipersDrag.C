/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
                                Copyright 2013-     Graz University of Technology
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "TangKuipersDrag.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(TangKuipersDrag, 0);

addToRunTimeSelectionTable
(
    forceModel,
    TangKuipersDrag,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
TangKuipersDrag::TangKuipersDrag
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    verbose_(false),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    densityFieldName_(propsDict_.lookup("densityFieldName")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    dPrim_(1.0),
    rhoParticle_(0.0),
    interpolation_(false),
    treatExplicit_(false),
    implDEM_(false),
    splitImplicitExplicit_(false),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
    scale_(1.),
    g_(9.81),
    interpolation_void_(false),
    interpolation_vel_(false)
{
	    
    if (propsDict_.found("verbose")) verbose_=true;
    if (propsDict_.found("treatExplicit")) treatExplicit_=true;
    if (propsDict_.found("interpolation")) interpolation_=true;
    if (propsDict_.found("interpolation_void")) 
    {
        interpolation_void_=true;
        Info << "Using interpolated voidage values." << endl;
    }
    if (propsDict_.found("interpolation_vel")) 
    {
        interpolation_vel_=true;
        Info << "Using interpolated velocity values." << endl;
    }
    if (propsDict_.found("splitImplicitExplicit"))
    {
        Info << "will split implicit / explicit force contributions." << endl;
        splitImplicitExplicit_ = true;
        if(!interpolation_) 
        Info << "WARNING: will only consider fluctuating particle velocity in implicit / explicit force split!" << endl;
    }
    if (propsDict_.found("implDEM"))
    {
        treatExplicit_=false;
        implDEM_=true;
        Info << "Using implicit DEM drag formulation." << endl;
    }
    
    particleCloud_.checkCG(true);

    //Set the reference length scale from particle information
    // get viscosity field
    #ifdef comp
        const volScalarField nufField = particleCloud_.turbulence().mu()/rho_;
    #else
        const volScalarField& nufField = particleCloud_.turbulence().nu();
    #endif
    scalar nuf(0);
    scalar rho(0);

    //Append the field names to be probed
    particleCloud_.probeM().initialize(typeName, "TangKuipersDrag.logDat");
    particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("Urel");        //other are debug
    particleCloud_.probeM().scalarFields_.append("Rep");          //other are debug
    particleCloud_.probeM().scalarFields_.append("beta");                 //other are debug
    particleCloud_.probeM().scalarFields_.append("voidfraction");       //other are debug
    particleCloud_.probeM().writeHeader();

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TangKuipersDrag::~TangKuipersDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void TangKuipersDrag::setForce() const
{
    scale_=particleCloud_.cg();
    Info << "TangKuipersDrag using scale from liggghts cg = " << scale_ << endl;

    // get viscosity field
    #ifdef comp
        const volScalarField nufField = particleCloud_.turbulence().mu()/rho_;
    #else
        const volScalarField& nufField = particleCloud_.turbulence().nu();
    #endif

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    label cellI=0;
    scalar beta(0);

    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
	scalar Vs(0);
	scalar localPhiP(0);
	scalar vCell(0);
	
    scalar F0=0.;
    scalar G0=0.;

	vector UfluidFluct(0,0,0);
    vector UsFluct(0,0,0);
    vector dragExplicit(0,0,0);

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);

    #include "setupProbeModel.H"

    for(int index = 0;index <  particleCloud_.numberOfParticles(); index++)
    {
        //if(mask[index][0])
        //{
            cellI = particleCloud_.cellIDs()[index][0];
            drag = vector(0,0,0);
            Ufluid= vector(0,0,0);
            
            if (cellI > -1) // particle Found
            {
                if(interpolation_)
                {
	                position     = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid       = UInterpolator_.interpolate(position,cellI);
                    //Ensure interpolated void fraction to be meaningful
                    if(voidfraction>1.00) voidfraction = 1.00;
                    if(voidfraction<0.10) voidfraction = 0.10;
                }
		else if(interpolation_void_)
                {
					position = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid = U_[cellI];//UInterpolator_.interpolate(position,cellI);
                    //Ensure interpolated void fraction to be meaningful
                    if(voidfraction>1.00) voidfraction = 1.00;
                    if(voidfraction<0.10) voidfraction = 0.10;
                }
		else if(interpolation_vel_)
                {
					position = particleCloud_.position(index);
                    voidfraction = particleCloud_.voidfraction(index);//voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid = UInterpolator_.interpolate(position,cellI);
                    //Ensure interpolated void fraction to be meaningful
                    if(voidfraction>1.00) voidfraction = 1.00;
                    if(voidfraction<0.10) voidfraction = 0.10;
                }
                else
                {
					voidfraction = voidfraction_[cellI];
                    Ufluid = U_[cellI];
                }
           
                Us = particleCloud_.velocity(index);
                Ur = Ufluid-Us;
                ds = 2*particleCloud_.radius(index);
                dPrim_ = ds/scale_;
                nuf = nufField[cellI];
                rho = rho_[cellI];

                magUr = mag(Ur);
				Rep = 0;
                Vs = ds*ds*ds*M_PI/6;
                localPhiP = 1-voidfraction+SMALL;
                vCell = U_.mesh().V()[cellI];
                
                 if (magUr > 0)
                {
                    // calc particle Re Nr
                    Rep = dPrim_*voidfraction*magUr/(nuf+SMALL);

                    // calc model coefficient F0
                    F0 = 10.f * localPhiP / voidfraction / voidfraction
                       + voidfraction * voidfraction * ( 1.0+1.5*Foam::sqrt(localPhiP) );

                    // calc model coefficient G0
					G0 =  Rep * ( 0.11 * localPhiP * ( 1.0 + localPhiP )
		    		- 0.00456 / voidfraction / voidfraction / voidfraction / voidfraction
		    		+ ( 0.169 * voidfraction + 0.0644 / voidfraction / voidfraction / voidfraction / voidfraction ) * Foam::pow(Rep,-0.343));


                    // calc model coefficient beta
                    beta =  18.*nuf*rho/(dPrim_*dPrim_)
                              *voidfraction * (F0 + G0);

                    
                    // calc particle's drag
                    drag = Vs * beta * Ur; //total drag force!

                    if (modelType_=="B")
                        drag /= voidfraction;
                        
                    //Split forces
                    if(splitImplicitExplicit_)
                    {
                        UfluidFluct  = Ufluid - U_[cellI];
                        UsFluct      = Us     - UsField_[cellI];
                        dragExplicit = Vs * beta * (UfluidFluct - UsFluct); //explicit part of force
    
                        if (modelType_=="B")
                            dragExplicit /= voidfraction;
                    }
                }

                //if( verbose_ ) //&& index>100 && index < 105)
               if(verbose_ && index < 2)  
               {
                    Pout << "index / cellI = " << index << tab << cellI << endl;
                    Pout << "position = " << particleCloud_.position(index) << endl;
                    Pout << "Us = " << Us << endl;
                    Pout << "Ur = " << Ur << endl;
                    Pout << "ds = " << ds << endl;
                    Pout << "rho = " << rho << endl;
                    Pout << "nuf = " << nuf << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "voidfraction_[cellI]: " << voidfraction_[cellI] << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "F0 / G0: " << F0 << tab << G0 << endl;
                    Pout << "beta:   " << beta << endl;
                    Pout << "drag = " << drag << endl;
                    
                    if(splitImplicitExplicit_)
                    {
                        Pout << "UfluidFluct = " << UfluidFluct << endl;
                        Pout << "UsFluct = " << UsFluct << endl;
                        Pout << "dragExplicit = " << dragExplicit << endl;
                    }

                    Pout << "\nTangKuipers drag model settings: treatExplicit " << treatExplicit_ << tab
                                    << "verbose: " << verbose_ << tab
                                    << "dPrim: " << dPrim_  << tab
                                    << "interpolation: " << interpolation_  << tab
                                    << endl << endl;
                }

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    vValues.append(drag);           //first entry must the be the force
                    vValues.append(Ur);                
                    sValues.append(Rep);
                    sValues.append(beta);
                    sValues.append(voidfraction);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }

            // set force on particle
            if(treatExplicit_) for(int j=0;j<3;j++) expForces()[index][j] += drag[j];
            else   //implicit treatment, taking explicit force contribution into account
            {
               for(int j=0;j<3;j++) 
               { 
                    impForces()[index][j] += drag[j] - dragExplicit[j]; //only consider implicit part!
                    expForces()[index][j] += dragExplicit[j];
               }
            }


            // set Cd
            if(implDEM_)
            {
                for(int j=0;j<3;j++) fluidVel()[index][j]=Ufluid[j];

                if (modelType_=="B")
                    Cds()[index][0] = Vs*beta/voidfraction;
                else
                    Cds()[index][0] = Vs*beta;

            }else{
                for(int j=0;j<3;j++) DEMForces()[index][j] += drag[j];
            }
        //}
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
