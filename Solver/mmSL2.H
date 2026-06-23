//==============================================================================
//  mmSL2.H — L2 error (numerical vs manufactured), uses mmsParams.H.
//==============================================================================
{
    scalar tNow  = runTime.value();
    scalar decay = Foam::exp(-tNow/tau_mms);

    scalar sumE = 0.0, sumL = 0.0, sumV = 0.0;
    forAll(mesh.C(), cellI)
    {
        scalar z   = mesh.C()[cellI].z() - z0_mms;
        scalar shp = Foam::cos(kz_mms*z);
        scalar TeEx = T0_mms + Ae_mms*shp*decay;
        scalar TlEx = T0_mms + Al_mms*shp*decay;
        scalar V = mesh.V()[cellI];
        sumE += Foam::sqr(Te[cellI] - TeEx) * V;
        sumL += Foam::sqr(Tl[cellI] - TlEx) * V;
        sumV += V;
    }
    reduce(sumE, sumOp<scalar>()); reduce(sumL, sumOp<scalar>()); reduce(sumV, sumOp<scalar>());
    scalar L2Te = Foam::sqrt(sumE/sumV);
    scalar L2Tl = Foam::sqrt(sumL/sumV);
    Info<< "MMS  t=" << tNow << "  L2(Te)=" << L2Te << "  L2(Tl)=" << L2Tl << endl;
    if (Pstream::master())
    {
        std::ofstream os("mms_error.dat", std::ios::app);
        os << tNow << " " << L2Te << " " << L2Tl << "\n";
    }
}