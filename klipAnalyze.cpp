



#include <mx/math/geo.hpp>

#include <mx/improc/fitsFile.hpp>

#include <mx/improc/imageFilters.hpp>
#include <mx/improc/eigenImage.hpp>
#include <mx/improc/eigenCube.hpp>
#include <mx/improc/ds9Interface.hpp>

#include <mx/gslInterpolation.hpp>
#include <mx/fileUtils.hpp>

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
 
      
      vec.push_back( mx::convertFromString<typeT>( str.substr(p0, p1-p0)));
      
      p0 = p1+1;
      //while(p0 < str.size() && isspace(str[p0])) ++p0;
      p1 = p0;
   }
   
   return vec;
}

template<typename floatT>
struct klipAnalyze
{
   std::vector<floatT> pas;
   std::vector<floatT> contrasts;
   std::vector<floatT> seps;
   std::vector<int> nmodes;
   floatT regminr;
   floatT regmaxr;
   floatT qthresh;
   int inclrefn;
   floatT mindpx;
   
   floatT fixedPA;
   floatT fixedSep;
   
   std::vector<floatT> stds;
   
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
      
   }
   
   void processHeader( mx::improc::fitsHeader & head )
   {
      if(pas.size() == 0)
      {
         pas = convertFromStringVector<floatT>(head["FAKEPA"].String());
      }
      if(contrasts.size() == 0)
      {
         contrasts = convertFromStringVector<floatT>(head["FAKECONT"].String());
      }
      if(seps.size() == 0)
      {
         seps = convertFromStringVector<floatT>(head["FAKESEP"].String());
      }
      if(nmodes.size() == 0)
      {
         nmodes = convertFromStringVector<int>(head["NMODES"].String());
      }
      if(regminr == -1)
      {
         regminr = mx::convertFromString<floatT>(head["REGMINR"].String());
      }
      if(regmaxr == -1)
      {
         regmaxr = mx::convertFromString<floatT>(head["REGMAXR"].String());
      }
      if(qthresh == -1)
      {
         qthresh = mx::convertFromString<floatT>(head["QTHRESH"].String());
      }
      if(inclrefn == -1)
      {
         inclrefn = head["INCLREFN"].Int();
      }
      if(mindpx == -1)
      {
         mindpx = head["MINDPX"].Value<floatT>();
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
    
   floatT getPAfromFileName(const std::string & fname)
   {
      int p1 = fname.size()-1;
      
      while(fname[p1] != '.') --p1;
      --p1;
      while(fname[p1] != '.') --p1;
      
      
      int p0 = p1-1;
      p1 += 3;

      while( isdigit(fname[p0]) ) --p0;
      ++p0;
      

      return mx::convertFromString<floatT>(fname.substr(p0, p1-p0));
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
   
      seps = {21.6};
     // cubeGaussUnsharpMask(ims, 20.0);
      cubeGaussSmooth(ims, 6.0);
      
      mx::improc::eigenImage<float> im, stdIm, mask;
   
      mask.resize(ims.rows(), ims.cols());
      mask.setConstant(1.0);
      
      floatT maskx = 0.5*(ims.cols()-1) - seps[0] * sin( mx::math::dtor(pas[0]) );
      floatT masky = 0.5*(ims.rows()-1) + seps[0] * cos( mx::math::dtor(pas[0]) );
      floatT maskr = 20;
      
      mx::improc::maskCircle(mask, maskx, masky, maskr);
            
      
      //ds9(mask);
      
      mx::improc::eigenCube<float> stdImc;
   
      mx::improc::stddevImageCube(stdImc, ims, mask, regminr, regmaxr, true);
   
      cubeGetMaxInMask(stds, ims, mask, 0);
      
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
      //positivePlanet();
   
     // cubeGaussUnsharpMask(ims, 20.0);
     // cubeGaussSmooth(ims, 6.0);
      mx::improc::eigenImage<float> im, stdIm, mask;
   
      mask.resize(ims.rows(), ims.cols());
      mask.setConstant(1.0);
      
      floatT maskx = 0.5*(ims.cols()-1) - seps[0] * sin( pas[0]*3.14159/180.0);
      floatT masky = 0.5*(ims.rows()-1) + seps[0] * cos( pas[0]*3.14159/180.0);
      
      //floatT maskx = 0.5*(ims.cols()-1) - fixedSep * sin( fixedPA*3.14159/180.0);
      //floatT masky = 0.5*(ims.rows()-1) + fixedSep * cos( fixedPA*3.14159/180.0);
      
      mx::improc::eigenImage<float> rIm, qIm;
      rIm.resize(ims.rows(), ims.cols());
      qIm.resize(ims.rows(), ims.cols());
      
      float xcen = 0.5*(ims.rows()-1);
      float ycen = 0.5*(ims.cols()-1);
      mx::improc::radAngImage( rIm, qIm, xcen, ycen );
      
      std::vector<size_t> maskdx =  mx::improc::annulusIndices( rIm, qIm, xcen, ycen, 16, 26, pas[0]-30+90, pas[0]+30+90);
                                                    
      
      floatT maskr = 20;
      
      //mx::improc::maskCircle(mask, maskx, masky, maskr);
            
      mx::improc::applyMask( mask, maskdx, 0);
      
      //ds9(mask);
      
      mx::improc::eigenCube<float> stdImc;
   
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
         std::cout << "0 0 0\n";
      }
   }
   
   void processFile( const std::string & fname, floatT sep = -1, bool parsePA = false )
   {
      initialize();
//       if(sep != -1) seps = {sep, sep};
//       if(parsePA)
//       {
//          floatT pa = getPAfromFileName(fname);
//          pas = {pa,pa};
//       }
      
      analyzeFilePP(fname);
      output(basename(fname.c_str()));
   
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

   std::vector<std::string> files = mx::getFileNames("/home/jrmales/findr/bpic0001", "output", "",".fits");
   
   //std::vector<std::string> files = mx::getFileNames("/home/jrmales/Data/Magellan/VisAO/2014.04.15/GQLup_ha_sdi_25s/bot/cent/findr/run1/", "output_97.1579", ".fits");
   
  
   
   
   #pragma omp parallel for
   for(int i=0; i<files.size(); ++i)
   {
      klipAnalyze<float> ka;
      
      #pragma omp critical
      std::cerr << basename(files[i].c_str()) << " " << i+1 << "/" << files.size() << "\n";
      ka.processFile(files[i]);//, 90.92, true);
   }
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


