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
#include "noEModel.H"
#include "addToRunTimeSelectionTable.H"

#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noEModel, 0);

addToRunTimeSelectionTable(eModel, noEModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
noEModel::noEModel
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    eModel(dict, sm)	
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noEModel::~noEModel()
{}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //
 
void noEModel::execute()
{}

const volScalarField noEModel::vField() const
{
	volScalarField vField_
	(
		IOobject
		(
		   "vField",
		   particleCloud_.mesh().time().timeName(),
		   particleCloud_.mesh(),
		   IOobject::NO_READ,
		   IOobject::NO_WRITE
		),	
		particleCloud_.mesh(),
		dimensionedScalar("zero", dimPower/dimCurrent, 0.0)
	);

	return vField_;    
}

void noEModel::postFlow()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
