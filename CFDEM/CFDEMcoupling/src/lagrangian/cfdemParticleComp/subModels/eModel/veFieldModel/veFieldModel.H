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
    Catalytic chemistry without mass transfer resistances

\*---------------------------------------------------------------------------*/

#ifndef veFieldModel_H
#define veFieldModel_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "eModel.H"

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class veFieldModel Declaration
\*---------------------------------------------------------------------------*/


class veFieldModel
:
    public eModel
{
private:

    dictionary propsDict_;

    //- Name of the voltage
    word voltageFieldName_; //CFD side
    const volScalarField& vField_;


    //- Name of the electric field
    word electricFieldName_; //CFD side
    const volVectorField& eField_;

    word particleVFieldName_; //DEM side
    word particleEFieldName_; //DEM side

    //- per particle fields
    mutable double **partV_;
    mutable double **partE_;
   
    // - Swicth to communicate fields to particles
    const bool usePartVE_;
    const bool interpolation_;


    virtual void giveData();

    // Private Member Functions
        // - Allocate arrays
        virtual void allocateMyArrays() const;

public:

    //- Runtime type information
    TypeName("veFieldModel");

    // Constructors

        //- Construct from components
        veFieldModel
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        ~veFieldModel();


    // Member Functions

        // - Calculate source terms
        virtual void execute();

	// - return electric potential for each particle
	virtual const volScalarField vField() const;

        virtual void postFlow();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
