/*---------------------------------------------------------------------------*\
License

    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this code.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "veFieldModel.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(veFieldModel, 0);

addToRunTimeSelectionTable(eModel, veFieldModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
veFieldModel::veFieldModel
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    eModel(dict, sm),
    propsDict_(dict.subDict(typeName + "Props")),

    // CFD side of fields
    voltageFieldName_(propsDict_.lookup("voltageFieldName")),
    vField_(sm.mesh().lookupObject<volScalarField> (voltageFieldName_)),
    electricFieldName_(propsDict_.lookup("electricFieldName")),
    eField_(sm.mesh().lookupObject<volVectorField> (electricFieldName_)),
    interpolation_(propsDict_.lookupOrDefault<bool>("interpolation",false)),

    // DEM - Particle fields
    partV_(NULL),
    partE_(NULL),
    usePartVE_(propsDict_.lookupOrDefault<bool>("useParticleVE",true)),
    particleVFieldName_(propsDict_.lookupOrDefault<word>("particleVFieldName", "pV")),
    particleEFieldName_(propsDict_.lookupOrDefault<word>("particleEFieldName", "pE"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

veFieldModel::~veFieldModel()
{
	//delete partV_;
	//delete partE_;
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void veFieldModel::allocateMyArrays() const
{
    particleCloud_.dataExchangeM().allocateArray(partV_, 0., 1);
    particleCloud_.dataExchangeM().allocateArray(partE_, 0., 3);
}

void veFieldModel::execute()
{
   //Info<<"execute electricModel"<<endl;
   allocateMyArrays();
   label cellI(0);
   scalar Vfluid(0.);
   vector Efluid(0,0,0);
	
   interpolationCellPoint<scalar> VInterpolator_(vField_);
   interpolationCellPoint<vector> EInterpolator_(eField_);

   // Assigning voltage and electric field from CFD to DEM - on each particle in the simulations 
   for(int index=0; index < particleCloud_.numberOfParticles(); ++index) 
   {   
	//Info<<"execute function in veEModel 2 "<<endl;
	cellI = particleCloud_.cellIDs()[index][0];
	if(cellI>=0)
	{
		if(interpolation_)
                {
                    vector position = particleCloud_.position(index);
                    Vfluid = VInterpolator_.interpolate(position,cellI);
                    Efluid = EInterpolator_.interpolate(position,cellI);
                }
                else
                {
                    Efluid = eField_[cellI];
                    Vfluid = vField_[cellI];
                }

		partV_[index][0]= Vfluid;
		partE_[index][0]=Efluid[0];
		partE_[index][1]=Efluid[1];
		partE_[index][2]=Efluid[2];       
		//Info<<"execute function in veEModel 3 :  "<<partV_[index][0]<<" , "<< partE_[index][0]<<" , "<<partE_[index][1]<< " ,"<< partE_[index][2]<< endl;
	}
   }  
}

const volScalarField veFieldModel::vField() const
{
	return vField_;
}


void veFieldModel::giveData()
{
    Info << "In the giveData function of eModel" << endl;
    particleCloud_.dataExchangeM().giveData(particleVFieldName_, "scalar-atom", partV_);
    particleCloud_.dataExchangeM().giveData(particleEFieldName_, "vector-atom", partE_);
}


void veFieldModel::postFlow()
{

  //Info<<"postFlow electricModel"<<endl;
  if(usePartVE_)
  {
    //Info<<"in the usePartVE_ if"<<endl;
   label cellI;
   scalar Vfluid(0.0);
   scalar Vpart(0.0);
   vector Efluid(0,0,0);
   
   interpolationCellPoint<scalar> VInterpolator_(vField_);
   interpolationCellPoint<vector> EInterpolator_(eField_);

   for( int index=0; index < particleCloud_.numberOfParticles(); ++index) 
   {   
	cellI = particleCloud_.cellIDs()[index][0];
	if(cellI>=0)
	{
		if(interpolation_)
                {
			vector position = particleCloud_.position(index);
                        Vfluid = VInterpolator_.interpolate(position,cellI);
                   	Efluid = EInterpolator_.interpolate(position,cellI);
                }
                else
                {
                       Vfluid = vField_[cellI];
                       Efluid = eField_[cellI];
                }

		partV_[index][0]= Vfluid;
		partE_[index][0]=Efluid[0];
		partE_[index][1]=Efluid[1];
		partE_[index][2]=Efluid[2];  	
		//Info<<"execute function in veEModel 3 :  "<<partV_[index][0]<<" , "<< partE_[index][0]<<" , "<<partE_[index][1]<< " ,"<< partE_[index][2]<< endl;

	}
   }
    //Info<<"Before giveData partV in postFlow veFieldModel"<<endl;
    giveData();   
    //Info<<"in the usePartVE_ end of if"<<endl;
  } 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
