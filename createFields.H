Info<< "Reading field PsiReal\n" << endl;

volScalarField PsiReal
(
    IOobject
    (
        "PsiReal",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field PsiImg\n" << endl;

volScalarField PsiImg
(
    IOobject
    (
        "PsiImg",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField absPsi
(
    IOobject
    (
        "absPsi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading diffusivity DT\n" << endl;

volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimViscosity, Zero)
);

Info<< "Reading diffusivity DThbar\n" << endl;

// m/s
dimensionSet dimhbar(1, 2, -1, 0, 0, 0, 0); // [kg m s K mol A Cd]

volScalarField DThbar
(
    IOobject
    (
        "DThbar",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimhbar, Zero)
);

Info<< "Reading mass\n" << endl;

volScalarField mass
(
    IOobject
    (
        "mass",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMass, Zero)
);


Info  << "dimhbar : " << dimhbar << endl;
if (!DThbar.headerOk())
{
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    Info  << "dimhbar : " << dimhbar << endl;
    DT = dimensionedScalar("DT", dimViscosity, transportProperties);
    DThbar = dimensionedScalar("DThbar", dimhbar , transportProperties);
    mass = dimensionedScalar("mass", dimMass , transportProperties);
    Info << "mass : " << mass << endl;
//    Info << "DThbar/(2*mass) : " << DThbar.value()/(2*mass.value()) << endl;
}

#include "createFvOptions.H"
