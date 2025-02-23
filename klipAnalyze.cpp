

#include <mx/math/geo.hpp>

#include <mx/ioutils/fits/fitsFile.hpp>
using namespace mx::fits;

#include <mx/improc/imageUtils.hpp>
#include <mx/improc/imageFilters.hpp>
#include <mx/improc/eigenImage.hpp>
#include <mx/improc/eigenCube.hpp>
using namespace mx::improc;

#include <mx/math/gslInterpolation.hpp>
#include <mx/ioutils/fileUtils.hpp>

#include <mx/math/fit/fitGaussian.hpp>
#include <mx/math/geo.hpp>
using namespace mx::math;

#include <mx/app/application.hpp>
using namespace mx::app;

#include "include/gridProcess.hpp"

using namespace mx::improc;

template <typename eigenImageT>
void imageGaussUnsharpMask( eigenImageT &im, typename eigenImageT::Scalar fwhm )
{
    eigenImageT fim;

    mx::improc::filterImage(
        fim, im, mx::improc::gaussKernel<eigenImageT, 2>( fwhm ) ); //, 0.5*(im.cols()-1) - fwhm*4);

    im = ( im - fim );
}

template <typename eigenCubeT>
void cubeGaussUnsharpMask( eigenCubeT &imc, typename eigenCubeT::Scalar fwhm )
{
#pragma omp parallel
    {
        typename eigenCubeT::imageT im;

#pragma omp for
        for( int i = 0; i < imc.planes(); ++i )
        {
            im = imc.image( i );
            imageGaussUnsharpMask( im, fwhm );
            imc.image( i ) = im;
        }
    }
}

template <typename eigenImageT>
void imageAzBoxUnsharpMask( eigenImageT &im, typename eigenImageT::Scalar radW, typename eigenImageT::Scalar azW )
{
    eigenImageT fim;

    filterImage( fim, im, azBoxKernel<eigenImageT>( radW, azW ) ); //, 0.5*(im.cols()-1) - 4*radW);

    im = ( im - fim );
}

template <typename eigenCubeT>
void cubeAzBoxUnsharpMask( eigenCubeT &imc, typename eigenCubeT::Scalar radW, typename eigenCubeT::Scalar azW )
{
#pragma omp parallel
    {
        typename eigenCubeT::imageT im;

#pragma omp for
        for( int i = 0; i < imc.planes(); ++i )
        {
            im = imc.image( i );
            imageAzBoxUnsharpMask( im, radW, azW );
            imc.image( i ) = im;
        }
    }
}

template <typename eigenImageT>
void imageGaussSmooth( eigenImageT &im, typename eigenImageT::Scalar fwhm )
{
    eigenImageT fim;

    mx::improc::filterImage(
        fim, im, mx::improc::gaussKernel<eigenImageT, 4>( fwhm ) ); //, 0.5*(im.cols()-1) - fwhm*4);

    im = fim;
}

template <typename eigenCubeT>
void cubeGaussSmooth( eigenCubeT &imc, typename eigenCubeT::Scalar fwhm )
{
#pragma omp parallel
    {
        typename eigenCubeT::imageT im;

#pragma omp for
        for( int i = 0; i < imc.planes(); ++i )
        {
            im = imc.image( i );
            imageGaussSmooth( im, fwhm );
            imc.image( i ) = im;
        }
    }
}

template <typename eigenImageT>
void imageAzBoxSmooth( eigenImageT &im, typename eigenImageT::Scalar radW, typename eigenImageT::Scalar azW )
{
    eigenImageT fim;

    filterImage( fim, im, azBoxKernel<eigenImageT>( radW, azW ) ); //, 0.5*(im.cols()-1) - 4*radW);

    im = fim;
}

template <typename eigenCubeT>
void cubeAzBoxSmooth( eigenCubeT &imc, typename eigenCubeT::Scalar radW, typename eigenCubeT::Scalar azW )
{
#pragma omp parallel
    {
        typename eigenCubeT::imageT im;

#pragma omp for
        for( int i = 0; i < imc.planes(); ++i )
        {
            im = imc.image( i );
            imageAzBoxSmooth( im, radW, azW );
            imc.image( i ) = im;
        }
    }
}

template <typename eigenImageT, typename eigenMaskT>
typename eigenImageT::Scalar imageGetMaxInMask( const eigenImageT &im, eigenMaskT &mask, int maskIncludeVal = 1 )
{
    typename eigenImageT::Scalar max = std::numeric_limits<typename eigenImageT::Scalar>::lowest();

    for( int j = 0; j < im.cols(); ++j )
    {
        for( int i = 0; i < im.rows(); ++i )
        {
            if( mask( i, j ) != maskIncludeVal )
                continue;

            if( im( i, j ) > max )
                max = im( i, j );
        }
    }

    return max;
}

template <typename eigenCubeT, typename eigenMaskT>
void cubeGetMaxInMask( std::vector<typename eigenCubeT::Scalar> &maxes,
                       eigenCubeT &imc,
                       eigenMaskT &mask,
                       int maskIncludeVal = 1 )
{
    maxes.resize( imc.planes() );

    for( int i = 0; i < imc.planes(); ++i )
    {
        maxes[i] = imageGetMaxInMask( imc.image( i ), mask, maskIncludeVal );
    }
}

template <typename eigenImageT, typename eigenMaskT>
typename eigenImageT::Scalar imageGetVarInMask( const eigenImageT &im, eigenMaskT &mask, int maskIncludeVal = 1 )
{

    std::vector<typename eigenImageT::Scalar> vals;

    for( int i = 0; i < im.rows(); ++i )
    {
        for( int j = 0; j < im.cols(); ++j )
        {
            if( mask( i, j ) != maskIncludeVal )
                continue;

            vals.push_back( im( i, j ) );
        }
    }

    if( vals.size() < 2 )
        return 0;

    return mx::math::vectorVariance( vals );
}

template <typename eigenCubeT, typename eigenMaskT>
void cubeGetVarInMask( std::vector<typename eigenCubeT::Scalar> &vars,
                       eigenCubeT &imc,
                       eigenMaskT &mask,
                       int maskIncludeVal = 1 )
{
    vars.resize( imc.planes() );

    for( int i = 0; i < imc.planes(); ++i )
    {
        vars[i] = imageGetVarInMask( imc.image( i ), mask, maskIncludeVal );
    }
}

template <typename realT>
int centroidImage( realT &x, realT &y, realT &A, realT &fwhm_x, realT &fwhm_y, realT &theta, eigenImage<realT> &im )
{
    realT Ag, xg, yg, xFWHM, yFWHM, angG;

    // Get a guess at the params
    mx::math::fit::guessGauss2D_ang<realT>( Ag, xg, yg, xFWHM, yFWHM, angG, im, 3, 15, 10, x, y );

    // Now fit it.
    mx::math::fit::fitGaussian2Dgen<realT> fg;
    fg.setGuess( 0, Ag, xg, yg, fwhm2sigma( xFWHM ), fwhm2sigma( yFWHM ), angG );

    fg.setArray( im.data(), im.rows(), im.cols() );
    fg.fit();

    x = fg.x0();
    y = fg.y0();
    A = fg.A();
    fwhm_x = sigma2fwhm( fg.sigma_x() );
    fwhm_y = sigma2fwhm( fg.sigma_y() );
    theta = fg.theta();

    return 0;
}

template <typename realT>
int centroidImageCube( std::vector<realT> &x,
                       std::vector<realT> &y,
                       std::vector<realT> &A,
                       std::vector<realT> &fwhm_x,
                       std::vector<realT> &fwhm_y,
                       std::vector<realT> &theta,
                       eigenCube<realT> &ims,
                       realT x0,
                       realT y0 )
{
    x.resize( ims.planes() );
    y.resize( ims.planes() );
    A.resize( ims.planes() );
    fwhm_x.resize( ims.planes() );
    fwhm_y.resize( ims.planes() );
    theta.resize( ims.planes() );

#pragma omp parallel
    {
        eigenImage<realT> im;
        realT tx, ty, tA, tfx, tfy, tt;

#pragma omp for
        for( int i = 0; i < ims.planes(); ++i )
        {
            tx = x0;
            ty = y0;

            im = ims.image( i );
            centroidImage( tx, ty, tA, tfx, tfy, tt, im );
            x[i] = tx;
            y[i] = ty;
            A[i] = tA;
            fwhm_x[i] = tfx;
            fwhm_y[i] = tfy;
            theta[i] = tt;
        }
    }

    return 0;
}

/// Convert a string to a std::vector of values
/** Parses the string based on the delimiter character template parameter.
 *
 * \param str is the std::string object to convert.
 *
 * \returns a vector of the converted numerical values.
 *
 * \tparam typeT the value type
 * \tparam delim is the delimiter character
 *
 */
template <typename typeT, char delim = ','>
inline std::vector<typeT> convertFromStringVector( const std::string &str )
{
    std::vector<typeT> vec;

    int p0, p1;

    p0 = 0;

    // while(p0 < str.size() && isspace(str[p0])) ++p0;
    p1 = p0;

    while( p1 < str.size() )
    {
        while( p1 < str.size() && str[p1] != delim )
            ++p1;

        vec.push_back( mx::ioutils::convertFromString<typeT>( str.substr( p0, p1 - p0 ) ) );

        p0 = p1 + 1;
        // while(p0 < str.size() && isspace(str[p0])) ++p0;
        p1 = p0;
    }

    return vec;
}

template <typename realT>
struct klipAnalyze : public application
{
    std::string m_file;
    int m_rows {0};
    int m_cols {0};
    std::vector<realT> m_planetX; ///< The x pixel coordinate of known planet signals
    std::vector<realT> m_planetY; ///< The y pixel coordinate of known planet signals
    std::vector<realT> m_planetR; ///< The radius of known planet signals (i.e. size to mask).  Default is 10.

    std::vector<realT> m_regminr; ///< The min. radii of the optimization regions.  Extracted from the header.
    std::vector<realT> m_regmaxr; ///< The min. radii of the optimization regions.  Extracted from the header.

    realT m_snrMinRad {0}; ///< The minimum radius for the SNR calculation.  Overrides min(m_regminr).
    realT m_snrMaxRad {0}; ///< The maximum radius for the SNR calculation.  Overrides max(m_regmaxr).

    realT m_hpfGaussFW {0};
    realT m_lpfGaussFW {0};

    std::vector<realT> m_pas;
    std::vector<realT> m_contrasts;
    std::vector<realT> m_seps;
    realT m_fakeRad {10};
    std::vector<int> m_nmodes;


    realT m_qthresh{ -1 };
    int m_inclrefn{ -1 };
    realT m_mindpx{ -1 };

    realT fixedPA;
    realT fixedSep;

    std::vector<realT> stds;
    std::vector<realT> As;
    std::vector<realT> drs;
    std::vector<realT> dqs;

    klipAnalyze()
    {
    }

    void initialize()
    {
        m_pas.clear();
        m_contrasts.clear();
        m_seps.clear();
        m_nmodes.clear();
        m_regminr = -1;
        m_regmaxr = -1;
        m_qthresh = -1;
        m_inclrefn = -1;
        m_mindpx = -1;

        stds.clear();
        As.clear();
        drs.clear();
        dqs.clear();
    }

    void setupConfig()
    {
        config.add( "file",
                    "f",
                    "file",
                    argType::Required,
                    "",
                    "file",
                    false,
                    "string",
                    "path to the input file" );

        config.add( "planet.X",
                    "X",
                    "planet.X",
                    argType::Required,
                    "planet",
                    "X",
                    false,
                    "vector<float>",
                    "x pixel coordinate of known planet signals" );
        config.add( "planet.Y",
                    "Y",
                    "planet.Y",
                    argType::Required,
                    "planet",
                    "Y",
                    false,
                    "vector<float>",
                    "y pixel coordinate of known planet signals" );
        config.add( "planet.R",
                    "R",
                    "planet.R",
                    argType::Required,
                    "planet",
                    "R",
                    false,
                    "vector<float>",
                    "The radius of known planet signals (i.e. size to mask)." );
        config.add( "snr.minRad",
                    "",
                    "snr.minRad",
                    argType::Required,
                    "snr",
                    "minRad",
                    false,
                    "float",
                    "inner radius in pixels of annulus for SNR calculation. Overrides optimization region from header." );
        config.add( "snr.maxRad",
                    "",
                    "snr.maxRad",
                    argType::Required,
                    "snr",
                    "maxRad",
                    false,
                    "float",
                    "outer radius in pixels of annulus for SNR calculation. Overrides optimization region from header." );

        config.add( "filter.hpfGaussFW",
                    "",
                    "filter.hpfGaussFW",
                    argType::Required,
                    "filter",
                    "hpfGaussFW",
                    false,
                    "float",
                    "FWHM of Gaussian kernel for high-pass filtering (unsharp mask)" );

        config.add( "filter.lpfGaussFW",
                    "",
                    "filter.lpfGaussFW",
                    argType::Required,
                    "filter",
                    "lpfGaussFW",
                    false,
                    "float",
                    "FWHM of Gaussian kernel for low-pass filtering (smoothing)" );
    }

    void loadConfig()
    {
        config(m_file, "file");
        config(m_planetX, "planet.X");
        config(m_planetY, "planet.Y");
        config(m_planetR, "planet.R");

        config(m_snrMinRad, "snr.minRad");

        config(m_snrMaxRad, "snr.maxRad");

        config(m_hpfGaussFW, "filter.hpfGaussFW");
        config(m_lpfGaussFW, "filter.lpfGaussFW");

    }

    void checkConfig()
    {
        if(m_planetX.size() > 0 || m_planetY.size() > 0 || m_planetR.size() > 0)
        {
            if(m_planetX.size() != m_planetY.size() || m_planetX.size() != m_planetR.size())
            {
                mxThrowException(mx::err::invalidconfig, "klipAnalyze::checkConfig", "planetX, planetY, and planetR must be same size");
            }
        }
    }

    void setup(int argc, char **argv)
    {
        return application::setup(argc, argv);
    }

    void processHeader( fitsHeader &head )
    {
        m_rows = head["NAXIS1"].Int();
        m_cols = head["NAXIS2"].Int();
        if( m_pas.size() == 0 && head.count("FAKEPA"))
        {
            m_pas = convertFromStringVector<realT>( head["FAKEPA"].String() );

            if( m_pas.size() == 0 )
            {
                std::cerr << "No fake PAs in header!\n";
                exit( -1 );
            }
        }
        if( m_contrasts.size() == 0 && head.count("FAKECONT"))
        {
            m_contrasts = convertFromStringVector<realT>( head["FAKECONT"].String() );

            if( m_contrasts.size() == 0 )
            {
                std::cerr << "No fake contrasts in header!\n";
                exit( -1 );
            }
        }
        if( m_seps.size() == 0 && head.count("FAKESEP"))
        {
            m_seps = convertFromStringVector<realT>( head["FAKESEP"].String() );

            if( m_seps.size() == 0 )
            {
                std::cerr << "No fake contrasts in header!\n";
                exit( -1 );
            }
        }

        for(size_t n = 0; n < m_seps.size(); ++n)
        {
            m_planetX.push_back(0.5*(m_rows-1.0)-m_seps[n]*sin(m_pas[n]*3.14159/180.));
            m_planetY.push_back(0.5*(m_cols-1.0)+m_seps[n]*cos(m_pas[n]*3.14159/180.));
            m_planetR.push_back(m_fakeRad);
        }


        if( m_nmodes.size() == 0 && head.count("NMODES"))
        {
            m_nmodes = convertFromStringVector<int>( head["NMODES"].String() );
        }

        if( m_regminr.size() == 0 && head.count("REGMINR"))
        {
            std::string rmr = head["REGMINR"].String();
            mx::ioutils::parseStringVector(m_regminr, rmr);
        }

        if( m_regmaxr.size() == 0 && head.count("REGMAXR"))
        {
            std::string rmr = head["REGMAXR"].String();
            mx::ioutils::parseStringVector(m_regmaxr, rmr);
        }

        if( m_qthresh == -1 && head.count("QTHRESH"))
        {
            m_qthresh = mx::ioutils::convertFromString<realT>( head["QTHRESH"].String() );
        }
        if( m_inclrefn == -1 && head.count("INCLREFN"))
        {
            m_inclrefn = head["INCLREFN"].Int();
        }
        if( m_mindpx == -1 && head.count("MINDPX"))
        {
            m_mindpx = head["MINDPX"].value<realT>();
        }
    }

    void positivePlanet()
    {
        for( int i = 0; i < m_contrasts.size(); ++i )
        {
            if( m_contrasts[i] < 0 )
            {
                m_contrasts.erase( m_contrasts.begin() + i );
                m_pas.erase( m_pas.begin() + i );
                m_seps.erase( m_seps.begin() + i );
                --i;
            }
        }
    }

    realT getPAfromFileName( const std::string &fname )
    {
        int p1 = fname.size() - 1;

        while( fname[p1] != '.' )
            --p1;
        --p1;
        while( fname[p1] != '.' )
            --p1;

        int p0 = p1 - 1;
        p1 += 3;

        while( isdigit( fname[p0] ) )
            --p0;
        ++p0;

        return mx::ioutils::convertFromString<realT>( fname.substr( p0, p1 - p0 ) );
    }

    void analyzeFilePP( const std::string &fname )
    {
        fitsFile<realT> ff;

        realT maskr = 15;
        realT smFWHM = 3.5;

        mx::improc::eigenCube<float> ims, proc;

        fitsHeader head;
        head.append( "FAKEPA" );
        head.append( "FAKECONT" );
        head.append( "FAKESEP" );
        head.append( "REGMINR" );
        head.append( "REGMAXR" );
        head.append( "NMODES" );
        head.append( "QTHRESH" );
        head.append( "MINDPX" );
        head.append( "INCLREFN" );

#pragma omp critical
        ff.read( ims, head, fname );

        // Extract the reduction params
        processHeader( head );

        // Remove any negative planets
        positivePlanet();

        // The center of the planet.
        realT cenx = 0.5 * ( ims.cols() - 1 ) - m_seps[0] * sin( dtor( m_pas[0] ) );
        realT ceny = 0.5 * ( ims.rows() - 1 ) + m_seps[0] * cos( dtor( m_pas[0] ) );

        std::cerr << "Found fake planet " << m_seps[0] << " " << m_pas[0] << " " << m_contrasts[0]
                  << " at pixel: " << cenx << " " << ceny << "\n";

        /*std::vector<realT> x;
        std::vector<realT> y;
        std::vector<realT> A;
        std::vector<realT> fwhm_x;
        std::vector<realT> fwhm_y;
        std::vector<realT> theta;*/

        // centroidImageCube( x, y, A, fwhm_x, fwhm_y, theta, ims, cenx, ceny);

        cubeGaussUnsharpMask( ims, 10.0 );
        cubeGaussSmooth( ims, smFWHM );

        mx::improc::eigenImage<float> im, stdIm, mask;

        mask.resize( ims.rows(), ims.cols() );
        mask.setConstant( 1.0 );

        mx::improc::maskCircle( mask, cenx, ceny, maskr, 0.0 );

        ff.write( "mask.fits", mask );

        mx::improc::eigenCube<float> stdImc;

        //mx::improc::stddevImageCube( stdImc, ims, mask, m_regminr, m_regmaxr, true );

        ff.write( "snrc.fits", stdImc );

        cubeGetMaxInMask( stds, stdImc, mask, 0 );

        /*
              //process results of fit
              As.resize(A.size());
              drs.resize(x.size());
              dqs.resize(y.size());

              realT msep, mq;
              realT tx, ty;

              for(int i=0;i< drs.size(); ++i)
              {
                 As[i] = A[i];

                 tx = x[i] - 0.5*(ims.cols()-1);
                 ty = y[i] - 0.5*(ims.rows()-1);
                 msep = sqrt(tx*tx + ty*ty);
                 drs[i] = msep - seps[0];

                 mq = angleMod( rtod(atan2( -tx, ty) ) );

                 dqs[i] = angleDiff(mq, pas[0]) ;
              }
              */
    }

#if 0
   void analyzeFileNP(const std::string & fname)
   {

      fitsFile<float> ff;

      mx::improc::eigenCube<float> ims, proc;

      fitsHeader head;
      head.append("FAKEPA");
      head.append("FAKECONT");
      head.append("FAKESEP");
      head.append("REGMINR");
      head.append("REGMAXR");
      head.append("NMODES");
      head.append("QTHRESH");
      head.append("MINDPX");
      head.append("INCLREFN");

    #pragma omp critical
      ff.read(ims, head, fname);

      processHeader(head);

      mx::improc::eigenImage<float> im, stdIm, mask;

      mask.resize(ims.rows(), ims.cols());
      mask.setConstant(1.0);

      //realT maskx = 0.5*(ims.cols()-1) - seps[0] * sin( pas[0]*3.14159/180.0);
      //realT masky = 0.5*(ims.rows()-1) + seps[0] * cos( pas[0]*3.14159/180.0);

      //realT maskx = 0.5*(ims.cols()-1) - fixedSep * sin( fixedPA*3.14159/180.0);
      //realT masky = 0.5*(ims.rows()-1) + fixedSep * cos( fixedPA*3.14159/180.0);

      mx::improc::eigenImage<float> rIm, qIm;
      rIm.resize(ims.rows(), ims.cols());
      qIm.resize(ims.rows(), ims.cols());

      float xcen = 0.5*(ims.rows()-1);
      float ycen = 0.5*(ims.cols()-1);
      mx::improc::radAngImage( rIm, qIm, xcen, ycen );

      std::vector<size_t> maskdx =  mx::improc::annulusIndices( rIm, qIm, xcen, ycen, 18, 30, pas[0]-10+90, pas[0]+10+90);


      //realT maskr = 20;

      //mx::improc::maskCircle(mask, maskx, masky, maskr);

      mx::improc::applyMask( mask, maskdx, 0);

      //ds9(mask);
      //exit(0);
      //mx::improc::eigenCube<float> stdImc;

     // mx::improc::stddevImageCube(stdImc, ims, mask, regminr, regmaxr, true);

      //cubeGetMaxInMask(stds, ims, mask, 0);


      cubeGetVarInMask(stds, ims, mask, 0);
   }
#endif

    void outputPP( bool maxOnly = false, const std::string &fname = "" )
    {
#pragma omp critical

        if( maxOnly )
        {
            realT max = std::numeric_limits<realT>::lowest();
            size_t maxn = 0;
            for( size_t n = 0; n < stds.size(); ++n )
            {
                if( stds[n] > max )
                {
                    maxn = n;
                    max = stds[n];
                }
            }

            std::cout << fname << " " << m_seps[0] << " " << m_pas[0] << " " << m_contrasts[0] << " " << m_qthresh
                      << " " << m_regminr << " " << m_regmaxr << " ";
            std::cout << m_mindpx << " " << m_inclrefn << " " << m_nmodes[maxn] << " " << stds[maxn] << "\n";
        }
        else
        {
            for( int i = 0; i < stds.size(); ++i )
            {
                std::cout << fname << " " << m_seps[0] << " " << m_pas[0] << " " << m_contrasts[0] << " " << m_qthresh
                          << " " << m_regminr << " " << m_regmaxr << " ";
                std::cout << m_mindpx << " " << m_inclrefn << " " << m_nmodes[i] << " " << stds[i] << "\n";
                // std::cout << As[i] << " " << drs[i] << " " << dqs[i] << "\n";
            }
        }
    }

#if 0
   void outputNP(const std::string & fname)
   {

      //realT min = 1e12;
      //int mini = -1;
    #pragma omp critical
      for(int i=0;i<stds.size(); ++i)
      {

         std::cout << seps[0] << " " << pas[0] << " " << contrasts[0] << " " <<  stds[i] << "\n";
      }

      //std::cerr << seps[0] << " " << pas[0] << " " << contrasts[0] << " " <<  stds[i] << "\n";
   }
#endif

#if 0
   void processFileNP( const std::string & fname, realT sep = -1, bool parsePA = false )
   {
      initialize();
//       if(sep != -1) seps = {sep, sep};
//       if(parsePA)
//       {
//          realT pa = getPAfromFileName(fname);
//          pas = {pa,pa};
//       }

      //analyzeFilePP(fname);
      analyzeFileNP(fname);
      outputNP(basename(fname.c_str()));

   }
#endif

    void processFilePP( const std::string &fname )
    {
        initialize();
        analyzeFilePP( fname );
        outputPP( false, basename( fname.c_str() ) );
    }
};

int main( int argc, char **argv )
{

    klipAnalyze<float> ka;

    ka.setup(argc, argv);

    float minRad = 50;
    float maxRad = 130;

    eigenCube<float> imc;
    fitsFile<float> ff;
    fitsHeader fh;

    ff.read( imc, fh, ka.m_file );


    ka.processHeader(fh);
    /*ka.m_snrMinRad=77;
    ka.m_snrMaxRad=107;*/

    //Get inner SNR annulus
    if(ka.m_snrMinRad <= 0)
    {
        ka.m_snrMinRad = *std::min_element(ka.m_regminr.begin(), ka.m_regminr.end());
    }

    //Get outer SNR annulus
    if(ka.m_snrMaxRad <= 0)
    {
        ka.m_snrMaxRad = *std::max_element(ka.m_regmaxr.begin(), ka.m_regmaxr.end());
    }

    eigenCube<float> maskCube;
    zeroNaNCube( imc, &maskCube );

    ff.write( "z.fits", imc );


    //cubeAzBoxUnsharpMask( imc, 3.0, 15 );

    if(ka.m_hpfGaussFW > 0)
    {
        cubeGaussUnsharpMask(imc, ka.m_hpfGaussFW);
    }

    if(ka.m_lpfGaussFW > 0)
    {
        cubeGaussSmooth( imc, ka.m_lpfGaussFW);
    }

    // cubeAzBoxSmooth(imc, (float) 1.5, (float) 2.75);

    eigenImage<float> planetMask;
    planetMask.resize( imc.rows(), imc.cols() );
    planetMask.setConstant( 1.0 );

    for(size_t n = 0; n < ka.m_planetX.size(); ++ n)
    {
        maskCircle( planetMask, ka.m_planetX[n], ka.m_planetY[n], ka.m_planetR[n], 0.0 );
    }

    ff.write( "planetMask.fits", planetMask );

    eigenCube<float> fullMaskCube;
    fullMaskCube.resize(maskCube.rows(), maskCube.cols(), maskCube.planes());

    for( int p = 0; p < maskCube.planes(); ++p )
    {
        maskCube.image( p) = -1*(maskCube.image(p) - 1);
        fullMaskCube.image( p ) = maskCube.image(p)* planetMask;
    }

    ff.write("analyzeMaskCube.fits", fullMaskCube);

    eigenCube<float> snrc;

    stddevImageCube( snrc, imc, fullMaskCube, ka.m_snrMinRad, ka.m_snrMaxRad, true );


    zeroNaNCube( snrc );

    for( int p = 0; p < maskCube.planes(); ++p )
    {
        snrc.image(p) *= maskCube.image(p);
    }

    ff.write( "snrc.fits", snrc );

    float best1sig = 1;

    eigenImage<float> mask;
    for( int p = 0; p < snrc.planes(); ++p )
    {
        // float sig1 = cont/(snrc.image(p)*(-1*(mask-1))).maxCoeff();
        // if(sig1 < best1sig) best1sig = sig1;
        mask = fullMaskCube.image(p);
        std::cerr << p << " ";
        if(ka.m_nmodes.size() == snrc.planes())
        {
            std::cerr << ka.m_nmodes[p] << " ";
        }
        std::cerr << ( snrc.image( p ) * ( -1 * ( mask - 1 ) ) ).maxCoeff() << '\n';
    }


}
