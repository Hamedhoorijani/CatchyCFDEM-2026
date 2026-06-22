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

    Description
    Electrification mode turned off

\*---------------------------------------------------------------------------*/

#ifndef noEModel_H
#define noEModel_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "eModel.H"

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noEModel Declaration
\*---------------------------------------------------------------------------*/


class noEModel
:
    public eModel
{
private:

   
public:

    //- Runtime type information
    TypeName("off");

    // Constructors

        //- Construct from components
        noEModel
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        ~noEModel();


    // Member Functions

        // - Calculate source terms
        virtual void execute();

	// - Return volScalarField of electric potential
	virtual const volScalarField vField() const;

        virtual void postFlow();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
