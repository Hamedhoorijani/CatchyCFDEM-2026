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
#include "eModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<eModel> eModel::New
(
    const dictionary& dict,
    cfdemCloudEnergy& sm,
    word eModelType
)
{
    Info<< "Selecting eModel "
         << eModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(eModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "eModel::New(const dictionary&, cfdemCloudEnergy&): "
            << endl
            << "    unknown eModelType type "
            << eModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid eModel types are :"
            << endl;

        Info<< dictionaryConstructorTablePtr_->toc()
            << abort(FatalError);
    }

    return autoPtr<eModel>(cstrIter()(dict,sm));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
