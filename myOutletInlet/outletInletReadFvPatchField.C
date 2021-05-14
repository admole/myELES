/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "volFields.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "outletInletReadFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

outletInletReadFvPatchScalarField::outletInletReadFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    phiName_("phi"),
    FieldName_("pIn"),
    GradFieldName_("pGradIn"),
    valueField_(p.size()),
    gradField_(p.size()),
    inlet_outlet_(0)
{
//    this->refValue() = *this;
//    this->refGrad() = Zero;
//    this->valueFraction() = 0.0;
}


outletInletReadFvPatchScalarField::outletInletReadFvPatchScalarField
(
    const outletInletReadFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    FieldName_(ptf.FieldName_),
    GradFieldName_(ptf.GradFieldName_),
    valueField_(ptf.valueField_, mapper),
    gradField_(ptf.gradField_, mapper),
    inlet_outlet_(ptf.inlet_outlet_)
{}


outletInletReadFvPatchScalarField::outletInletReadFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    FieldName_(dict.lookupOrDefault<word>("FieldName", "pIn")),
    GradFieldName_(dict.lookupOrDefault<word>("GradFieldName", "pGradIn")),
    valueField_(p.size()),
    gradField_(p.size()),
    inlet_outlet_(dict.lookupOrDefault("inletOutlet", true))
{
    //patchType() = dict.lookupOrDefault<word>("patchType", word::null);

    const fvMesh& mesh_(patch().boundaryMesh().mesh());

    const volScalarField valuein
            (
                    IOobject
                            (
                                    FieldName_,
                                    mesh_.time().timeName(),
                                    mesh_,
                                    IOobject::MUST_READ,
                                    IOobject::AUTO_WRITE,
                                    true // registery
                            ),
                    mesh_
            );

    const volScalarField gradin
            (
                    IOobject
                            (
                                    GradFieldName_,
                                    mesh_.time().timeName(),
                                    mesh_,
                                    IOobject::MUST_READ,
                                    IOobject::AUTO_WRITE,
                                    true // registery
                            ),
                    mesh_
            );

    const label patchId = mesh_.boundaryMesh().findPatchID(this->patch().name());

    valueField_ = valuein.boundaryField()[patchId]; //->internalField();
    gradField_ = gradin.boundaryField()[patchId].snGrad(); //->internalField();

    refValue() = Field<scalar>("outletValue", dict, p.size());

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
                (
                        scalarField("value", dict, p.size())
                );
    }
    else
    {
        fvPatchScalarField::operator=(refValue());
    }

    refGrad() = Zero;
    valueFraction() = 0.0;
}


outletInletReadFvPatchScalarField::outletInletReadFvPatchScalarField
(
    const outletInletReadFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    phiName_(ptf.phiName_),
    FieldName_(ptf.FieldName_),
    GradFieldName_(ptf.GradFieldName_),
    valueField_(ptf.valueField_),
    gradField_(ptf.gradField_),
    inlet_outlet_(ptf.inlet_outlet_)
{}


outletInletReadFvPatchScalarField::outletInletReadFvPatchScalarField
(
    const outletInletReadFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    phiName_(ptf.phiName_),
    FieldName_(ptf.FieldName_),
    GradFieldName_(ptf.GradFieldName_),
    valueField_(ptf.valueField_),
    gradField_(ptf.gradField_),
    inlet_outlet_(ptf.inlet_outlet_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void outletInletReadFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh &mesh = patch().boundaryMesh().mesh();
    const label patchId = mesh.boundaryMesh().findPatchID(this->patch().name());
    const volScalarField &valuein = db().objectRegistry::lookupObject<volScalarField>(FieldName_);
    const volScalarField &gradin = db().objectRegistry::lookupObject<volScalarField>(GradFieldName_);

    valueField_ = valuein.boundaryField()[patchId]; //->internalField();
    gradField_ = gradin.boundaryField()[patchId].snGrad(); //->internalField();


    Info<< "outletInlet Reading from RANS "
        << endl;

    const Field<scalar>& phip =
        patch().template lookupPatchField<surfaceScalarField, scalar>
        (
            phiName_
        );

    if (inlet_outlet_)
    {
        valueFraction() = pos0(phip);
    }
    else
    {
        valueFraction() = 1.0;
    }
    refGrad() = gradField_;
    refValue() = valueField_;


    mixedFvPatchField<scalar>::updateCoeffs();
}


void outletInletReadFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    refValue().writeEntry("outletValue", os);
    writeEntry("value", os);
    os.writeKeyword("inletOutlet") << inlet_outlet_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        makePatchTypeField
        (
                fvPatchScalarField,
                outletInletReadFvPatchScalarField
        );

// ************************************************************************* //


} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
