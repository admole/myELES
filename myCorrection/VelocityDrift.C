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

#include "VelocityDrift.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrix.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(VelocityDrift, 0);
    addToRunTimeSelectionTable
    (
        option,
        VelocityDrift,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::VelocityDrift::drift(fvMatrix<vector>& eqn)
{

    // access to velocity equation terms
    const scalarField& vol = mesh_.V();
    const volVectorField& U = eqn.psi();
    // scalarField& diag = eqn.diag();
    vectorField& source = eqn.source();

    //const dimensionedScalar& readNu = db().lookupObject<IOdictionary>("transportProperties").lookup("nu");

    //const volVectorField U = mesh_.lookupObject<volVectorField>("U");
    const volVectorField& U_les = mesh_.lookupObject<volVectorField>("U_LES");
    const volScalarField& k_ = mesh_.lookupObject<volScalarField>("turbulenceProperties:k");
    const volScalarField& eps_ = mesh_.lookupObject<volScalarField>("turbulenceProperties:epsilon");

    const dimensionedScalar smallv("smallv", dimensionSet(0,1,-1,0,0,0,0), 1e-06);
    const volVectorField flowDir_ = U_les / max(mag(U_les), smallv);         // should this be U or U_les (or combination)
    const volVectorField position_ = mesh_.C();

    const vector b_size = b_max - b_min;
    const scalar max_dist = mag(b_size);
    const scalar dt_ = mesh_.time().deltaTValue();

    const scalar t = mesh_.time().value();
    const scalar ramp = min((t - rampStart)/(rampEnd - rampStart), 1.0);


    forAll(cells_, i)
    {
        label celli = cells_[i];

        vector dist_normal(0.0, 0.0, 0.0);
        vector dists(0.0, 0.0, 0.0);

        for(int j=0; j<3; j++)
        {
            if (flowDir_[celli][j] > 0.0)
            {
                dist_normal[j] = b_max[j] - position_[celli][j];
                dists[j] = dist_normal[j]/flowDir_[celli][j];
            }
            else if (flowDir_[celli][j] < 0.0)
            {
                dist_normal[j] = position_[celli][j] - b_min[j];
                dists[j] = dist_normal[j]/flowDir_[celli][j];
            }
            else if ((flowDir_[celli][0] == 0.0) && (flowDir_[celli][1] == 0.0) && (flowDir_[celli][2] == 0.0))
            {
                dist_normal[j] = min(b_max[j] - position_[celli][j], position_[celli][j] - b_min[j]);
                dists[j] = dist_normal[j];
            }
            else
            {
                dists[j] = max_dist;
            }

            // dists[j] = dist_normal[j]/max(flowDir_[celli][j], SMALL);
        }

        scalar dist = cmptMin(cmptMag(dists));
        scalar norm_dist = dist / max_dist;
        scalar wi = 1 - norm_dist;
        scalar timescale = max(C1_*k_[celli]/eps_[celli], dt_);
        // scalar timescale = 1.0;
        // scalar timescale = dt_;
        vector drift = ramp * wi * (U[celli] - U_les[celli]) / (timescale) * vol[celli];
        source[celli] += drift;
        
    }

    Info<< type() << " " << name_ << " corrected U to LES Mean " << endl;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::VelocityDrift::VelocityDrift
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

void Foam::fv::VelocityDrift::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    drift(eqn);
}

void Foam::fv::VelocityDrift::correct
        (
                volVectorField& U
        )
{
    U.correctBoundaryConditions();
}

void Foam::fv::VelocityDrift::writeData(Ostream& os) const
{
    dict_.writeEntry(name_, os);
}


bool Foam::fv::VelocityDrift::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.readEntry("bMin", b_min);
        coeffs_.readEntry("bMax", b_max);
        coeffs_.readEntry("rampStart", rampStart);
        coeffs_.readEntry("rampEnd", rampEnd);
        coeffs_.readEntry("C1", C1_);           //TODO: Read if present???

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
