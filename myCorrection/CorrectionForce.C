/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "CorrectionForce.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrix.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(CorrectionForce, 0);
    addToRunTimeSelectionTable
    (
        option,
        CorrectionForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::CorrectionForce::addDamping(fvMatrix<vector>& eqn)
{
    // Note: we want to add
    //      deltaU/deltaT
    // where deltaT is a local time scale:
    //  U/(cbrt of volume)
    // Since directly manipulating the diagonal we multiply by volume.


    // access to velocity equation terms
    // const scalarField& vol = mesh_.V();
    const volVectorField& U = eqn.psi();
    // scalarField& diag = eqn.diag();
    vectorField& source = eqn.source();

    //const dimensionedScalar& readNu = db().lookupObject<IOdictionary>("transportProperties").lookup("nu");

    //const volVectorField U = mesh_.lookupObject<volVectorField>("U");
    const volVectorField& U_les = mesh_.lookupObject<volVectorField>("U_LES");

    const dimensionedScalar smallv("smallv", dimensionSet(0,1,-1,0,0,0,0), 1e-06);
    const volVectorField flowDir_ = U_les / max(mag(U_les), smallv);         // should this be U or U_les (or combination)
    const volVectorField position_ = mesh_.C();

    //const vector b_min(10.0,0.0,-1.57);       // TODO: replace with reading from dictionary
    //const vector b_max(12.0,2.0,1.57);
    const vector b_size = b_max - b_min;
    const scalar max_dist = std::max(std::max(b_size.x(), b_size.y()), b_size.z()); // better way of finding max of vector


    forAll(cells_, i)
    {
        label celli = cells_[i];

        vector dist_normal(0.0, 0.0, 0.0);
        vector dists(0.0, 0.0, 0.0);

        for(int j=0; j<3; j++)
        {
            if (flowDir_[celli][j] > 0.0){ dist_normal[j] = b_max[j] - position_[celli][j];}
            else                         { dist_normal[j] = position_[celli][j] - b_min[j];}

                dists[j] = dist_normal[j]/max(flowDir_[celli][j], SMALL);
        }

        scalar dist = std::min(std::min(dists.x(), dists.y()), dists.z()); // better way of finding min of vector
        scalar wi = 1 - dist / max(max_dist, SMALL);
        // scalar relax2 = max(0.02*k/eps, dt)

        source[celli] += wi * (U_les[celli] - U[celli]) / 1000; // divide by relaxation term
    }

    // TODO: Do I need to update the boundary conditions here?

    Info<< type() << " " << name_ << " corrected U to LES Mean " << endl;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::CorrectionForce::CorrectionForce
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::CorrectionForce::constrain
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    addDamping(eqn);
}


void Foam::fv::CorrectionForce::writeData(Ostream& os) const
{
    dict_.writeEntry(name_, os);
}


bool Foam::fv::CorrectionForce::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.readEntry("bMin", b_min);
        coeffs_.readEntry("bMax", b_max);

        if (!coeffs_.readIfPresent("UNames", fieldNames_))
        {
            fieldNames_.resize(1);
            fieldNames_.first() = coeffs_.lookupOrDefault<word>("U", "U");
        }

        applied_.setSize(fieldNames_.size(), false);

        return true;
    }

    return false;
}


// ************************************************************************* //
