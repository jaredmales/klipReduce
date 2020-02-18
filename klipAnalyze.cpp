



#include <mx/math/geo.hpp>

#include <mx/improc/fitsFile.hpp>

#include <mx/improc/imageFilters.hpp>
#include <mx/improc/eigenImage.hpp>
#include <mx/improc/eigenCube.hpp>
#include <mx/improc/ds9Interface.hpp>

#include <mx/gslInterpolation.hpp>
#include <mx/ioutils/fileUtils.hpp>

#include <mx/math/fit/fitGaussian.hpp>
#include <mx/math/geo.hpp>
using namespace mx::math;

#include "include/gridProcess.hpp"

using namespace mx::improc;

template<typename eigenImageT>
void imageGaussUnsharpMask(eigenImageT & im, typename eigenImageT::Scalar fwhm)
{
   eigenImageT fim;

   mx::improc::filterImage(fim, im, mx::improc::gaussKernel<eigenImageT,2>(fwhm));//, 0.5*(im.cols()-1) - fwhm*4);
         
   im = (im-fim);
}

template<typename eigenCubeT>
void cubeGaussUnsharpMask(eigenCubeT & imc, typename eigenCubeT::Scalar fwhm)
{
   #pragma omp parallel
   {
      typename eigenCubeT::imageT im;

      #pragma omp for
      for(int i=0; i< imc.planes(); ++i)
      {
         im = imc.image(i);
         imageGaussUnsharpMask(im, fwhm);
         imc.image(i) = im;
      }
   }
}

template<typename eigenImageT>
void imageAzBoxUnsharpMask(eigenImageT & im, typename eigenImageT::Scalar radW, typename eigenImageT::Scalar azW)
{
   eigenImageT fim;

   filterImage(fim, im, azBoxKernel<eigenImageT>( radW, azW));//, 0.5*(im.cols()-1) - 4*radW);
   
   im = (im-fim);
}

template<typename eigenCubeT>
void cubeAzBoxUnsharpMask(eigenCubeT & imc, typename eigenCubeT::Scalar radW, typename eigenCubeT::Scalar azW)
{
   #pragma omp parallel
   {
      typename eigenCubeT::imageT im;

      #pragma omp for
      for(int i=0; i< imc.planes(); ++i)
      {
         im = imc.image(i);
         imageAzBoxUnsharpMask(im, radW, azW);
         imc.image(i) = im;
      }
   }
}

template<typename eigenImageT>
void imageGaussSmooth(eigenImageT & im, typename eigenImageT::Scalar fwhm)
{
   eigenImageT fim;

   mx::improc::filterImage(fim, im, mx::improc::gaussKernel<eigenImageT,4>(fwhm));//, 0.5*(im.cols()-1) - fwhm*4);
         
   im = fim;
}

template<typename eigenCubeT>
void cubeGaussSmooth(eigenCubeT & imc, typename eigenCubeT::Scalar fwhm)
{
   #pragma omp parallel
   {
      typename eigenCubeT::imageT im;
      
      #pragma omp for
      for(int i=0; i< imc.planes(); ++i)
      {
         im = imc.image(i);
         imageGaussSmooth(im, fwhm);
         imc.image(i) = im;
      }
   }
}

template<typename eigenImageT>
void imageAzBoxSmooth(eigenImageT & im, typename eigenImageT::Scalar radW, typename eigenImageT::Scalar azW)
{
   eigenImageT fim;

   filterImage(fim, im, azBoxKernel<eigenImageT>( radW, azW));//, 0.5*(im.cols()-1) - 4*radW);
   
   im = fim;
}

template<typename eigenCubeT>
void cubeAzBoxSmooth(eigenCubeT & imc, typename eigenCubeT::Scalar radW, typename eigenCubeT::Scalar azW)
{
   #pragma omp parallel
   {
      typename eigenCubeT::imageT im;

      #pragma omp for
      for(int i=0; i< imc.planes(); ++i)
      {
         im = imc.image(i);
         imageAzBoxSmooth(im, radW, azW);
         imc.image(i) = im;
      }
   }
}


template<typename eigenImageT, typename eigenMaskT>
typename eigenImageT::Scalar imageGetMaxInMask( const eigenImageT & im,
                                                eigenMaskT & mask,
                                                int maskIncludeVal = 1 )
{
   typename eigenImageT::Scalar max = 0;
   
   for(int i=0; i<im.rows(); ++i)
   {
      for(int j=0; j<im.cols(); ++j)
      {
         if( mask(i,j) != maskIncludeVal ) continue;
         
         if( im(i,j) > max) max = im(i,j);
      }
   }
   
   return max;
}

template<typename eigenCubeT, typename eigenMaskT>
void cubeGetMaxInMask( std::vector<typename eigenCubeT::Scalar> & maxes,
                       eigenCubeT & imc,
                       eigenMaskT & mask,
                       int maskIncludeVal = 1 )
{
   maxes.resize( imc.planes());
   
   for(int i=0; i<imc.planes(); ++i)
   {
      maxes[i] = imageGetMaxInMask( imc.image(i), mask, maskIncludeVal);
   }
   
}

template<typename eigenImageT, typename eigenMaskT>
typename eigenImageT::Scalar imageGetVarInMask( const eigenImageT & im,
                                                   eigenMaskT & mask,
                                                   int maskIncludeVal = 1 )
{
   
   std::vector<typename eigenImageT::Scalar> vals;
   
   for(int i=0; i<im.rows(); ++i)
   {
      for(int j=0; j<im.cols(); ++j)
      {
         if( mask(i,j) != maskIncludeVal ) continue;
   
         vals.push_back(im(i,j));
      }
   }

   if(vals.size() < 2) return 0;
   
   return mx::math::vectorVariance(vals);
}

template<typename eigenCubeT, typename eigenMaskT>
void cubeGetVarInMask( std::vector<typename eigenCubeT::Scalar> & vars,
                       eigenCubeT & imc,
                       eigenMaskT & mask,
                       int maskIncludeVal = 1 )
{
   vars.resize( imc.planes());
   
   for(int i=0; i<imc.planes(); ++i)
   {
      vars[i] = imageGetVarInMask( imc.image(i), mask, maskIncludeVal);
   }
   
}

template<typename realT>
int centroidImage( realT & x,
                   realT & y,
                   realT & A,
                   realT & fwhm_x,
                   realT & fwhm_y,
                   realT & theta,
                   eigenImage<realT> & im
                 )
{
   realT Ag, xg, yg, xFWHM, yFWHM, angG;
   
   //Get a guess at the params
   mx::math::fit::guessGauss2D_ang<realT>( Ag, xg, yg, xFWHM, yFWHM, angG, im, 3, 15, 10, x, y);

   //Now fit it.
   mx::math::fit::fitGaussian2Dgen<realT> fg;
   fg.setGuess(0, Ag, xg, yg, fwhm2sigma(xFWHM), fwhm2sigma(yFWHM), angG);
   
   fg.setArray(im.data(), im.rows(), im.cols());
   fg.fit();
   
   x = fg.x0();
   y = fg.y0();
   A = fg.A();
   fwhm_x = sigma2fwhm(fg.sigma_x());
   fwhm_y = sigma2fwhm(fg.sigma_y());
   theta = fg.theta();
   
   return 0;
}

template<typename realT>
int centroidImageCube( std::vector<realT> & x,
                       std::vector<realT> & y,
                       std::vector<realT> & A,
                       std::vector<realT> & fwhm_x,
                       std::vector<realT> & fwhm_y,
                       std::vector<realT> & theta,
                       eigenCube<realT> & ims,
                       realT x0,
                       realT y0
                     )
{
   x.resize(ims.planes());
   y.resize(ims.planes());
   A.resize(ims.planes());
   fwhm_x.resize(ims.planes());
   fwhm_y.resize(ims.planes());
   theta.resize(ims.planes());
   
   #pragma omp parallel
   {
      eigenImage<realT> im;
      realT tx, ty, tA, tfx, tfy, tt;

      #pragma omp for
      for(int i=0; i< ims.planes(); ++i)
      {
         tx = x0;
         ty = y0;
         
         im = ims.image(i);
         centroidImage(tx, ty, tA, tfx, tfy, tt, im);
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
template<typename typeT, char delim=','>  inline
std::vector<typeT> convertFromStringVector(const std::string & str)
{
   std::vector<typeT> vec;
   
   int p0, p1;
   
   p0 = 0;
   
   //while(p0 < str.size() && isspace(str[p0])) ++p0;
   p1 = p0;
   
   while(p1 < str.size())
   {
      while(p1 < str.size() && str[p1] != delim) ++p1;
 
      
      vec.push_back( mx::ioutils::convertFromString<typeT>( str.substr(p0, p1-p0)));
      
      p0 = p1+1;
      //while(p0 < str.size() && isspace(str[p0])) ++p0;
      p1 = p0;
   }
   
   return vec;
}

template<typename realT>
struct klipAnalyze
{
   std::vector<realT> pas;
   std::vector<realT> contrasts;
   std::vector<realT> seps;
   std::vector<int> nmodes;
   realT regminr;
   realT regmaxr;
   realT qthresh;
   int inclrefn;
   realT mindpx;
   
   realT fixedPA;
   realT fixedSep;
   
   std::vector<realT> stds;
   std::vector<realT> As;
   std::vector<realT> drs;
   std::vector<realT> dqs;
   
   
   mx::improc::ds9Interface ds9;
   
   klipAnalyze()
   {
      regminr = -1;
      regmaxr = -1;
      qthresh = -1;
      inclrefn = -1;
      mindpx = -1;
   }
   
   void initialize()
   {
      pas.clear();
      contrasts.clear();
      seps.clear();
      nmodes.clear();
      regminr = -1;
      regmaxr = -1;
      qthresh = -1;
      inclrefn = -1;
      mindpx = -1;
      
      stds.clear();
      As.clear();
      drs.clear();
      dqs.clear();
      
   }
   
   void processHeader( mx::improc::fitsHeader & head )
   {
      if(pas.size() == 0)
      {
         pas = convertFromStringVector<realT>(head["FAKEPA"].String());
         
         if(pas.size() == 0)
         {
            std::cerr << "No fake PAs in header!\n";
            exit(-1);
         }

         
      }
      if(contrasts.size() == 0)
      {
         contrasts = convertFromStringVector<realT>(head["FAKECONT"].String());
         
         if(contrasts.size() == 0)
         {
            std::cerr << "No fake contrasts in header!\n";
            exit(-1);
         }
      }
      if(seps.size() == 0)
      {
         seps = convertFromStringVector<realT>(head["FAKESEP"].String());
         
         if(seps.size() == 0)
         {
            std::cerr << "No fake contrasts in header!\n";
            exit(-1);
         }

      }
      if(nmodes.size() == 0)
      {
         nmodes = convertFromStringVector<int>(head["NMODES"].String());
      }
      if(regminr == -1)
      {
         regminr = mx::ioutils::convertFromString<realT>(head["REGMINR"].String());
      }
      if(regmaxr == -1)
      {
         regmaxr = mx::ioutils::convertFromString<realT>(head["REGMAXR"].String());
      }
      if(qthresh == -1)
      {
         qthresh = mx::ioutils::convertFromString<realT>(head["QTHRESH"].String());
      }
      if(inclrefn == -1)
      {
         inclrefn = head["INCLREFN"].Int();
      }
      if(mindpx == -1)
      {
         mindpx = head["MINDPX"].Value<realT>();
      }
   }
   
   void positivePlanet()
   {
      for(int i=0; i< contrasts.size(); ++i)
      {
         if(contrasts[i] < 0)
         {
            contrasts.erase( contrasts.begin()+i);
            pas.erase( pas.begin() + i);
            seps.erase( seps.begin() + i);
            --i;
         }
      }
   }
    
   realT getPAfromFileName(const std::string & fname)
   {
      int p1 = fname.size()-1;
      
      while(fname[p1] != '.') --p1;
      --p1;
      while(fname[p1] != '.') --p1;
      
      
      int p0 = p1-1;
      p1 += 3;

      while( isdigit(fname[p0]) ) --p0;
      ++p0;
      

      return mx::ioutils::convertFromString<realT>(fname.substr(p0, p1-p0));
   }
   
      
   void analyzeFilePP(const std::string & fname)
   {
      mx::improc::fitsFile<float> ff;

      mx::improc::eigenCube<float> ims, proc;
   
      mx::improc::fitsHeader head;
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
      
      positivePlanet();

      
      
      realT cenx = 0.5*(ims.cols()-1) - seps[0] * sin( dtor(pas[0]) );
      realT ceny = 0.5*(ims.rows()-1) + seps[0] * cos( dtor(pas[0]) );
      

      std::vector<realT> x;
      std::vector<realT> y;
      std::vector<realT> A;
      std::vector<realT> fwhm_x;
      std::vector<realT> fwhm_y;
      std::vector<realT> theta;


      centroidImageCube( x, y, A, fwhm_x, fwhm_y, theta, ims, cenx, ceny);
   

     // cubeGaussUnsharpMask(ims, 20.0);
      cubeGaussSmooth(ims, 6.0);
      
            
      mx::improc::eigenImage<float> im, stdIm, mask;
   
      mask.resize(ims.rows(), ims.cols());
      mask.setConstant(1.0);
            
      
      realT maskr = 20;
      
      mx::improc::maskCircle(mask, cenx, ceny, maskr);
            
      
      //ds9(mask);
      
      mx::improc::eigenCube<float> stdImc;
   
      mx::improc::stddevImageCube(stdImc, ims, mask, regminr, regmaxr, true);
   
      cubeGetMaxInMask(stds, ims, mask, 0);
      

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
      
   }
   
   
   void analyzeFileNP(const std::string & fname)
   {
      
      mx::improc::fitsFile<float> ff;

      mx::improc::eigenCube<float> ims, proc;
   
      mx::improc::fitsHeader head;
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
   
   
   void output(const std::string & fname)
   {
            
      #pragma omp critical
      for(int i=0;i<stds.size(); ++i)
      {
         std::cout << fname << " " << seps[0] << " " << pas[0] << " " << contrasts[0] << " " << qthresh << " " << regminr << " " << regmaxr << " ";
         std::cout << mindpx << " " << inclrefn << " " << nmodes[i] << " " << stds[i] << " ";
         std::cout << As[i] << " " << drs[i] << " " << dqs[i] << "\n";
      }
   }
   
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
   
   void processFile( const std::string & fname, realT sep = -1, bool parsePA = false )
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
   
   void processFilePP( const std::string & fname)
   {
      initialize();
      analyzeFilePP(fname);
      output(basename(fname.c_str()));
   
   }
   
};


int main()
{

   /*klipAnalyze<float> ka;
   ka.pas = {120};
   ka.seps = {46};
   ka.contrasts = {1e-5};
   
   ka.analyzeFilePP("/home/jrmales/Downloads/finim0000.fits");
   ka.output("");*/
   
   eigenCube<float> imc;
   fitsFile<float> ff;
   //ff.read(imc,"/home/jrmales/tmp/klip/contrast_0.00005/fake_28px_120deg_0000.fits");
   ff.read(imc, "finim0012.fits");
   
   //eigenImage<float> im50 = imc.image(9);
   
   for(int ii=0;ii<imc.rows(); ++ii)
   {
      for(int jj=0; jj<imc.cols(); ++jj)
      {
         for(int pp=0;pp<imc.planes(); ++pp)
         {
            if( !isnormal(imc.image(pp)(ii,jj)))  imc.image(pp)(ii,jj)=0;
         }
      }
   }
   
   //cubeGaussUnsharpMask(imc, 10.0);
   //cubeAzBoxUnsharpMask(imc, 3,15);
   cubeGaussUnsharpMask(imc, 10.0);
   cubeGaussSmooth(imc, (float) 4.8);
   //cubeAzBoxSmooth(imc, (float) 4, (float) 2);
   
   ds9Interface ds9(imc, 1);
   
   /*eigenImage<float> mask;
   mask.resize(imc.rows(), imc.cols());
   mask.setConstant(1.0);
   
   maskCircle( mask, 87, 97, 20 );
   
   eigenCube<float> snrc;
   stddevImageCube( snrc, imc, mask, 10, 80, true); 

   ds9(imc,2);
   ds9(mask, 3);

   ds9(snrc,4);*/
   
#if 0
   
   std::vector<std::string> files = mx::getFileNames("/home/jrmales/Data/Magellan/Clio/clio_20141202_03/bpic/findr/bpic39002/reduced", "output", "",".fits");

   
   #pragma omp parallel for
   for(int i=0; i<files.size(); ++i)
   {
      klipAnalyze<float> ka;
      
      #pragma omp critical
      std::cerr << basename(files[i].c_str()) << " " << i+1 << "/" << files.size() << "\n";
      ka.processFilePP(files[i]);//, 90.92, true);
   }
#endif 
#if 0

   std::vector<std::string> files = mx::getFileNames("/home/jrmales/Data/Magellan/Clio/clio_20141202_03/bpic39/neg2/", "finim", "",".fits");

   mx::improc::eigenCube<float> imc;
   mx::improc::fitsFile<float> ff;
   
   //ff.read(imc, files);
   
   //ff.write("neg.fits", imc);
   
   //return 0;
   
   //#pragma omp parallel for
   for(int i=0; i<files.size(); ++i)
   {
      klipAnalyze<float> ka;
      
      #pragma omp critical
      std::cerr << basename(files[i].c_str()) << " " << i+1 << "/" << files.size() << "\n";
      
      ka.processFile(files[i]);//, 90.92, true);
   }

#endif 
#if 0

   gridFile<double>("out", "/home/jrmales/Data/Magellan/Clio/clio_20141202_03/bpic/findr/bpic39002/reduced/results_nonan.txt");
   
#endif
   
   
//    //float pa = ka.getPAfromFileName(fname);
//    //ka.pas = {pa, pa};
//    ka.analyzeFile("/home/jrmales/Data/Magellan/VisAO/2014.04.15/GQLup_ha_sdi_25s/bot/cent/findr/run1/output_327.1585.fits");
//    
//    
// //    for(int i=0; i< ims.planes(); ++i)
// //    {
// //       im = ims.image(i);
// //       mx::stddevImage(stdIm, im, mask, 70,120, true);
// //       //imageGaussSmooth(im, 5);
// //       ims.image(i) = im; //im/stdIm*mask;
// //    }
// 
//    
//    
//    //for(int i=0; i<ka.stds.size(); ++i) std::cout << ka.stds[i] << "\n";
//    ka.output();
}


