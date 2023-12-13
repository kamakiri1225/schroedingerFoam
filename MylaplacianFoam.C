/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Application
    laplacianFoam

Group
    grpBasicSolvers

Description
    Laplace equation solver for a scalar quantity.

    \heading Solver details
    The solver is applicable to, e.g. for thermal diffusion in a solid.  The
    equation is given by:

    \f[
        \ddt{T}  = \div \left( D_T \grad T \right)
    \f]

    Where:
    \vartable
        T     | Scalar field which is solved for, e.g. temperature
        D_T   | Diffusion coefficient
    \endvartable

    \heading Required fields
    \plaintable
        T     | Scalar field which is solved for, e.g. temperature
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Laplace equation solver for a scalar quantity."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            // solve(fvm::ddt(T) - fvm::laplacian(DT, T));

            fvScalarMatrix PsiRealEqn
            (
                fvm::ddt(PsiImg) - fvm::laplacian(DThbar/(2*mass), PsiReal)
             //==
             //   fvOptions(PsiImg)
            );

            //add
            /*
                fvcにして陽解法にしている
                空間微分にクランクニコルソン法を適用する場合は、fvc::laplacian(DT, Ti)ではなく
                0.5*( fvm::laplacian(DT, Ti) + fvc::laplacian(DT, Ti) )
                のように、TnとToの平均を計算するようにする
            */
            fvScalarMatrix PsiImgEqn
            (
                fvm::ddt(PsiReal) + fvm::laplacian(DThbar/(2*mass), PsiImg)
             //==
             //   fvOptions(PsiReal)
            );

            //fvOptions.constrain(PsiRealEqn);
            PsiRealEqn.solve();
            //fvOptions.correct(PsiReal);

            //add 
            //fvOptions.constrain(PsiImgEqn);
            PsiImgEqn.solve();
            //fvOptions.correct(PsiImg);
            //add end
        }

        forAll(PsiReal,i)
        {
            absPsi[i] = sqr(pow(PsiReal[i],2) + pow(PsiImg[i],2) );
        }
        //Info << "absT : " << absT << endl;
        
        //#include "write.H"
        runTime.write(); //add
        runTime.printExecutionTime(Info);
    }
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
