/*---------------------------------------------------------------------------*\
License
    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Author: Alex Skillen.       alex.skillen@stfc.ac.uk

Reference. Accuracy and Efficiency Improvements in Synthetic Eddy Methods. A. Skillen A. Revell T. Craft IJHFF (2016) DOI: 10.1016/j.ijheatfluidflow.2016.09.008
\*---------------------------------------------------------------------------*/

#include "SEMBase.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMeshEntries.H"
#include "SEMspot.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTemplateTypeNameAndDebugWithName(IOList<symmTensor>, "symmTensorList", 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


SEMBase::SEMBase
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(p, iF),
    UFieldName_(iF.name()),
    UGradFieldName_(iF.name()),
    EpsFieldName_(iF.name()),
    RFieldName_(iF.name()),
    ranGen_(label(0)),
    meanField_(p.size()),
    UGradIn_(p.size()),
    curTimeIndex_(-1),
    RIn_(p.size()),
    epsIn_(p.size()),
    sigma_(p.size()),
    maxDelta_(p.size()),
    maxSigma_(0.1),
    embedded_(0),
    inlet_(1),
    outlet_(0),
    phiName_("phi")
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}



SEMBase::SEMBase
(
    const SEMBase& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<vector>(ptf, p, iF, mapper),
    UFieldName_(iF.name()),
    UGradFieldName_(iF.name()),
    EpsFieldName_(iF.name()),
    RFieldName_(iF.name()),
    ranGen_(label(0)),
    meanField_(ptf.meanField_, mapper),
    UGradIn_(ptf.UGradIn_, mapper),
    curTimeIndex_(-1),
    semBox_(ptf.semBox_),
    RIn_(ptf.RIn_, mapper),
    epsIn_(ptf.epsIn_, mapper),
    sigma_(ptf.sigma_, mapper),
    maxDelta_(ptf.maxDelta_, mapper),
    maxSigma_(ptf.maxSigma_),
    embedded_(ptf.embedded_),
    inlet_(ptf.inlet_),
    outlet_(ptf.outlet_),
    avgWindow_(ptf.avgWindow_),
    phiName_(ptf.phiName_)
{
}



SEMBase::SEMBase
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<vector>(p, iF),
    UFieldName_(dict.lookup("UFieldName")),
    UGradFieldName_(dict.lookup("UGradFieldName")),
    EpsFieldName_(dict.lookup("EpsFieldName")),
    RFieldName_(dict.lookup("RFieldName")),
    ranGen_(label(0)),
    meanField_(p.size()),
    UGradIn_(p.size()),
    //meanField_(this->db().lookupObject<volVectorField>("UIn")),
    curTimeIndex_(-1),
    semBox_(p.Cf()), 
    RIn_(p.size()),
    //RIn_(this->db().lookupObject<volSymmTensorField>("RIn")),
    epsIn_(p.size()),
    //epsIn_(this->db().lookupObject<volScalarField>("epsIn")),
    sigma_(p.size()),
    embedded_(dict.lookupOrDefault("embedded", true)),
    inlet_(dict.lookupOrDefault("inlet", true)),
    outlet_(dict.lookupOrDefault("outlet", false)),
    maxDelta_(p.size()),
    maxSigma_(dict.lookupOrDefault<scalar>("maxSigma", GREAT)),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))

{

    this->patchType() = dict.lookupOrDefault<word>("patchType", word::null);

    this->refValue() = Field<vector>("inletValue", dict, p.size());

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
                (
                        Field<vector>("value", dict, p.size())
                );
    }
    else
    {
        fvPatchField<vector>::operator=(this->refValue());
    }

    this->refGrad() = Zero;
    this->valueFraction() = 0.0;


//    meanField_ = Field<vector>("UIn", dict, p.size());
//    RIn_ = Field<symmTensor>("RIn", dict, p.size());
//    epsIn_ = Field<scalar>("epsIn", dict, p.size());
    //sigma_ = Field<vector>("sigma", dict, p.size());

    const fvMesh& mesh_(this->patch().boundaryMesh().mesh());

    const volVectorField Uin
            (
                    IOobject
                            (
                                    UFieldName_,
                                    mesh_.time().timeName(),
                                    mesh_,
                                    IOobject::MUST_READ,
                                    IOobject::AUTO_WRITE,
                                    true // registery
                            ),
                    mesh_
            );

    const volVectorField UGradin
            (
                    IOobject
                            (
                                    UGradFieldName_,
                                    mesh_.time().timeName(),
                                    mesh_,
                                    IOobject::MUST_READ,
                                    IOobject::AUTO_WRITE,
                                    true // registery
                            ),
                    mesh_
            );

    const volSymmTensorField Rin
            (
                    IOobject
                            (
                                    RFieldName_,
                                    mesh_.time().timeName(),
                                    mesh_,
                                    IOobject::MUST_READ,
                                    IOobject::AUTO_WRITE,
                                    true // registery
                            ),
                    mesh_
            );

    const volScalarField Epsin
            (
                    IOobject
                            (
                                    EpsFieldName_,
                                    mesh_.time().timeName(),
                                    mesh_,
                                    IOobject::MUST_READ,
                                    IOobject::AUTO_WRITE,
                                    true // registery
                            ),
                    mesh_
            );


    const label patchId = mesh_.boundaryMesh().findPatchID(this->patch().name());

    meanField_ = Uin.boundaryField()[patchId]; //->internalField();
    UGradIn_ = UGradin.boundaryField()[patchId].snGrad(); //->internalField();
    epsIn_ = Epsin.boundaryField()[patchId]; //->internalField();
    RIn_ = Rin.boundaryField()[patchId]; //->internalField();


}



SEMBase::SEMBase
(
    const SEMBase& ptf
)
:
    mixedFvPatchField<vector>(ptf),
    UFieldName_(ptf.UFieldName_),
    UGradFieldName_(ptf.UFieldName_),
    EpsFieldName_(ptf.EpsFieldName_),
    RFieldName_(ptf.RFieldName_),
    ranGen_(ptf.ranGen_),
    meanField_(ptf.meanField_),
    UGradIn_(ptf.UGradIn_),
    curTimeIndex_(-1),
    semBox_(ptf.semBox_), 
    RIn_(ptf.RIn_),
    epsIn_(ptf.epsIn_),
    sigma_(ptf.sigma_),
    embedded_(ptf.embedded_),
    inlet_(ptf.inlet_),
    outlet_(ptf.outlet_),
    maxDelta_(ptf.maxDelta_),
    maxSigma_(ptf.maxSigma_),
    spot_(ptf.spot_), 
    avgWindow_(ptf.avgWindow_),
    phiName_(ptf.phiName_)
{
}



SEMBase::SEMBase
(
    const SEMBase& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(ptf, iF),
    UFieldName_(ptf.UFieldName_),
    UGradFieldName_(ptf.UFieldName_),
    EpsFieldName_(ptf.EpsFieldName_),
    RFieldName_(ptf.RFieldName_),
    ranGen_(ptf.ranGen_),
    meanField_(ptf.meanField_),
    UGradIn_(ptf.UGradIn_),
    curTimeIndex_(ptf.curTimeIndex_),
    semBox_(ptf.semBox_),
    UBulk_(ptf.UBulk_),
    RIn_(ptf.RIn_),
    epsIn_(ptf.epsIn_),
    sigma_(ptf.sigma_),
    embedded_(ptf.embedded_),
    inlet_(ptf.inlet_),
    outlet_(ptf.outlet_),
    maxDelta_(ptf.maxDelta_),
    maxSigma_(ptf.maxSigma_),
    avgWindow_(ptf.avgWindow_),
    phiName_(ptf.phiName_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void SEMBase::initilise()
{

// FvPatchField<vector>::operator==(vector(0.0, 0.0, 0.0));

    if( embedded_ )
    {

        Info<< "Embedded SEM "
            << endl;

        if( this->db().time().value() != 0.0 )
        {
            //const fvMesh& mesh = dimensionedInternalField().mesh();
            const fvMesh &mesh = patch().boundaryMesh().mesh();
            const label patchId = mesh.boundaryMesh().findPatchID(this->patch().name());
            const volVectorField &Uin = db().objectRegistry::lookupObject<volVectorField>(UFieldName_);
            const volVectorField &UGradin = db().objectRegistry::lookupObject<volVectorField>(UGradFieldName_);
            const volScalarField &Epsin = db().objectRegistry::lookupObject<volScalarField>(EpsFieldName_);
            const volSymmTensorField &Rin = db().objectRegistry::lookupObject<volSymmTensorField>(RFieldName_);

            meanField_ = Uin.boundaryField()[patchId]; //->internalField();
            UGradIn_ = UGradin.boundaryField()[patchId].snGrad(); //->internalField();
            epsIn_ = Epsin.boundaryField()[patchId]; //->internalField();
            RIn_ = Rin.boundaryField()[patchId]; //->internalField();

            Info<< "SEM Reading from RANS "
                << endl;
        }
        else
        {
           //
        }
    }
    else
    {
        //
    }

    Info<< "compute bulk velocity "
        << endl;
    //compute bulk velocity of patch
    scalar patchTotArea( 0.0 );
    UBulk_ = 0.0;

    forAll(*this, facei)
    {
        UBulk_ += ( meanField_[facei] & this->patch().nf()()[facei] ) * this->patch().magSf()[facei];
        patchTotArea += this->patch().magSf()[facei];
    }

    reduce( UBulk_, sumOp<scalar>() );
    reduce( patchTotArea, sumOp<scalar>() );

    UBulk_ /= max( patchTotArea, SMALL);

    //get max dx dy dz of each cell adjoining the face
    //This is used for restricting minimum allowable eddy size.
    const fvMesh& mesh(this->patch().boundaryMesh().mesh());

    forAll (*this, facei)
    {
        const edgeList& cellEdges( mesh.cells()[this->patch().faceCells()[facei]].edges( mesh.faces() ) );
        scalarList cellEdgeLengths( cellEdges.size(), 0.0 );

        forAll (cellEdges, edgei)
        {
            cellEdgeLengths[edgei] = mag( cellEdges[edgei].vec( mesh.points() ) );
        }

        maxDelta_[facei] = max(cellEdgeLengths);
    }

    // set length-scales
    Info<< "Set Length-Scales "
        << endl;
    //scalar C=(isA<turbulentInletDFSEMFvPatchField>(*this)) ?  1.825742 : 1.0;
    scalar C = 1.0;

    forAll(*this, facei)
    {
        sigma_[facei] =
                vector
                        (
                                pow(mag(RIn_[facei].xx()), 1.5),
                                pow(mag(RIn_[facei].yy()), 1.5),
                                pow(mag(RIn_[facei].zz()), 1.5)
                        ) *C / epsIn_[facei];

        for( int i=0; i<3; i++ )
        {
            sigma_[facei][i] = max( sigma_[facei][i],  maxDelta_[facei] );
            sigma_[facei][i] = min( sigma_[facei][i],  maxSigma_ );
        }
    }

    Info<< "Set SEM box "
        << endl;
    // set SEM box size
    forAll(*this, facei)
    {
        vector inflate( sigma_[facei] - ( this->patch().nf()()[facei] & sigma_[facei] ) * this->patch().nf()()[facei] );

        for( int k=0; k<3; k++ )
        {
            semBox_.min()[k] = min( semBox_.min()[k], this->patch().Cf()[facei][k] - inflate[k] );
            semBox_.max()[k] = max( semBox_.max()[k], this->patch().Cf()[facei][k] + inflate[k] );
        }
    }
    
    reduce( semBox_.min(), minOp<vector>() );
    reduce( semBox_.max(), maxOp<vector>() );

    Info<< "Set Window "
        << endl;

    //set averaging window size
    avgWindow_ = cmptMax( max(sigma_) ) / max(mag(UBulk_), SMALL) * 5.0;

    reduce( avgWindow_, maxOp<scalar>() );

    Info<< "Window = "
            << avgWindow_
            << endl;
}

int SEMBase::numEddies()
{
    int numSpots = 0;
    const int maxSpots=100000;

    scalar minLength = GREAT;
    
    forAll(*this, facei)
    {
        minLength = min( minLength, cmptMin(sigma_[facei]));
    }

    reduce( minLength, minOp<scalar>() );

    scalar patchTotArea( gSum(this->patch().magSf()) );
    numSpots = 5.0 * patchTotArea / (2.0/3.0*3.14 * pow(minLength,2) + SMALL);

    reduce( numSpots, maxOp<scalar>() );
    label requested = numSpots;
    numSpots = min( numSpots, maxSpots );
    if( numSpots == maxSpots )
    {
        WarningIn("SEMBase::numEddies()")
            << "Maximum number of eddies reached\n"
            << "\tRequested: "
            << requested
            << endl
            << "\tLimited to: "
            << maxSpots 
            << endl;
    }
    else
    {
        Info<< "Number of eddies used: " 
            << numSpots
            << endl;
    }

    return numSpots; 
}



//---
void SEMBase::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchField<vector>::autoMap(m);
    meanField_.autoMap(m);
    UGradIn_.autoMap(m);
    RIn_.autoMap(m);
    epsIn_.autoMap(m);
    //sigma_.autoMap(m);
}


//---
void SEMBase::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<vector>::rmap(ptf, addr);

    const SEMBase& tiptf =
        refCast<const SEMBase >(ptf);

    meanField_.rmap(tiptf.meanField_, addr);
    UGradIn_.rmap(tiptf.UGradIn_, addr);
    RIn_.rmap(tiptf.RIn_, addr);
    epsIn_.rmap(tiptf.epsIn_, addr);
    //sigma_.rmap(tiptf.sigma_, addr);
}



//---
void SEMBase::advectPoints()
{
    //advect random points
    scalar dt = this->db().time().deltaTValue();

    for( int i=0; i<spot_.size(); i++ ) 
    {
        spot_[i]->origin() += spot_[i]->u()*dt;

        bool regen =  this->db().time().value() > spot_[i]->residenceTime();
        reduce( regen, orOp<bool>() );

        if(regen)
        {
            Info<< "Regenerating spot as time > "
                << spot_[i]->residenceTime()
                << endl;
            scalar origSize = mag( spot_[i]->sigma() );
            Info<< "Old spot size = "
                << origSize
                << endl;
            do
            {
                Info<< "Attempting regeneration"<< endl;
                spot_[i]->initialise(true);
            } while( mag( spot_[i]->sigma() ) > 1.1*origSize || mag( spot_[i]->sigma() ) < 0.9*origSize );
            Info<< "New spot size = "
                << mag( spot_[i]->sigma() )
                << endl;
        }
    }
}





void SEMBase::correctMass() 
{
    //correct bulk velocity of patch
    scalar patchTotArea = 0.0;
    scalar Uc = 0.0;

    forAll(*this, facei)
    {
        Uc += ( (*this)[facei] & this->patch().nf()()[facei] ) * this->patch().magSf()[facei];
        patchTotArea += this->patch().magSf()[facei];
    }

    reduce( Uc, sumOp<scalar>() );
    reduce( patchTotArea, sumOp<scalar>() );


    Uc /= max( patchTotArea, SMALL);

    forAll(*this, facei)
    { 
        (*this)[facei] *= mag(UBulk_)/max(mag(Uc), SMALL);
    } 
}



//---
void SEMBase::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    
    if (curTimeIndex_ != this->db().time().timeIndex() && this->db().time().value() > this->db().time().deltaTValue())
    {
        Info<< "SEM inlet / outlet patch = "
            << this->patch().name()
            << endl;

        this->advectPoints();

        this->updateU();

        this->correctMass();

        this->refGrad() = UGradIn_;

        label patchIndex = this->patch().index();
        //const label patchId = mesh.boundaryMesh().findPatchID(this->patch().name());
        const surfaceScalarField & phi = this->db().objectRegistry::lookupObject<surfaceScalarField>(phiName_);
        const scalarField & phip = phi.boundaryField()[patchIndex];

        if (inlet_ && outlet_)
        {
            //this->valueFraction() = 1.0 - pos0(phip);
            this->valueFraction() = 1.0 - pos0(meanField_ & this->patch().nf());
            Info<< "value fraction = "
                << this->valueFraction()
                << endl;

        }
        else if (inlet_)
        {
            this->valueFraction() = 1.0;
        }
        else if (outlet_)
        {
            this->valueFraction() = 0.0;
        }
        else
        {
            Info<< "WARNING: Both inlet and outlet switched off - defaulting to inlet/outlet"
                << endl;
            this->valueFraction() = 1.0 - pos0(phip);
        }

        curTimeIndex_ = this->db().time().timeIndex();
    }

    mixedFvPatchField<vector>::updateCoeffs();

}






void SEMBase::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("embedded") << embedded_ << token::END_STATEMENT << nl;
    os.writeKeyword("inlet") << inlet_ << token::END_STATEMENT << nl;
    os.writeKeyword("outlet") << outlet_ << token::END_STATEMENT << nl;
    os.writeKeyword("maxSigma") << maxSigma_ << token::END_STATEMENT << nl;
    os.writeKeyword("UFieldName") << UFieldName_ << token::END_STATEMENT << nl;
    os.writeKeyword("UGradFieldName") << UGradFieldName_ << token::END_STATEMENT << nl;
    os.writeKeyword("EpsFieldName") << EpsFieldName_ << token::END_STATEMENT << nl;
    os.writeKeyword("RFieldName") << RFieldName_ << token::END_STATEMENT << nl;
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    this->refValue().writeEntry("inletValue", os);

    this->writeEntry("value", os);

    //meanField_.writeEntry("UIn", os);
    //RIn_.writeEntry("RIn", os);
    //epsIn_.writeEntry("epsIn", os);
    //sigma_.writeEntry("sigma", os);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

    // template<class Type>

    void SEMBase::operator=
            (
                    const fvPatchField<vector>& ptf
            )
    {
        fvPatchField<vector>::operator=
                (
                        this->valueFraction()*this->refValue()
                        + (1 - this->valueFraction())*ptf
                );
    }


// ************************************************************************* //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
