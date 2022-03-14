// Solve the Momentum equation
tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(voidfraction*rho, U) + fvm::div(phi, U)
  + particleCloud.divVoidfractionTau(U, voidfraction)
  + fvm::Sp(Ksl, U)
 ==
    voidfraction*rho*g
  + fvOptions(voidfraction, rho, U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor() && (modelType=="B" || modelType=="Bfull"))
{
    UEqn == -fvc::grad(p) + Ksl*Us;
    UEqn.solve();
    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
else if (pimple.momentumPredictor())
{
    UEqn == -voidfraction*fvc::grad(p) + Ksl*Us;
    UEqn.solve();
    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
    Info<< "U max/min: " << gMax(U) << "/" << gMin(U) << endl;
