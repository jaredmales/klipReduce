
#include <iostream>



#define MX_APP_DEFAULT_configPathGlobal_env "KLIPREDUCE_GLOBAL_CONFIG"
#define MX_APP_DEFAULT_configPathLocal "klipReduce.conf"

#include <mx/gnuPlot.hpp>

#include <mx/app/application.hpp>


#include <mx/improc/ADIDerotator.hpp>
#include <mx/improc/KLIPreduction.hpp>

#include <libgen.h>

/// A program to run the mxlib KLIP pipeline
/**
  * 
  */
template<typename _realT, typename _evCalcT=double>
class klipReduce : public mx::application
{
public:
   typedef _realT realT;
   typedef _evCalcT evCalcT;
   
protected:

   //Basic Setup [HCIobservation]:
   std::string directory;
   std::string prefix;
   std::string extension;
   std::string fileList;
   
   //File Reading [HCIobservation]: 
   int deleteFront;
   int deleteBack;
   std::string qualityFile;
   realT qualityThreshold;
   bool thresholdOnly;
   
   
   std::string MJDKeyword;
   bool MJDisISO8601;
   int imSize;

   //Rotation setup [ADIobservation]
   std::string angleKeyword;
   float angleScale;
   float angleConstant;
   
   //Fake Planet Injection [ADIobservation]
   int fakeMethod;
   std::string fakeFileName; ///<FITS file containing the fake planet PSF to inject
   std::string fakeScaleFileName; ///< One-column text file containing a scale factor for each point in time.
   
   std::vector<realT> fakeSep; ///< Separation(s) of the fake planet(s)
   std::vector<realT> fakePA; ///< Position angles(s) of the fake planet(s)
   std::vector<realT> fakeContrast; ///< Contrast(s) of the fake planet(s)
     
   //Co-adding
   std::string coaddMethod;
   int coaddMaxImno;
   realT coaddMaxTime;
   
   //Masking
   std::string maskFile;
   
   //Pre-processing
   bool preProcess_beforeCoadd;
   bool preProcess_mask;
   bool preProcess_subradprof;
   realT preProcess_azUSM_azW;
   realT preProcess_azUSM_radW;
   realT preProcess_gaussUSM_fwhm;
   std::string preProcess_outputPrefix;
   bool preProcess_only;
   bool skipPreProcess;
   
   //KLIP parameters
   realT  minDPx;
   std::string excludeMethod;
   int includeRefNum;
   std::vector<int> Nmodes;
   std::vector<realT> minRadius;
   std::vector<realT> maxRadius;

   
   //Combination and Output
   
   bool noDerotate;
   
   std::string combineMethod;
   std::string weightFile;
   _realT sigmaThreshold;
   _realT minGoodFract {0.0};
   
   std::string outputFile;
   bool exactFName;
   std::string outputDir;
   
   //Output the individual PSF subtracted images
   bool outputPSFSub;
   std::string psfSubPrefix;
   
   //************************************//
   // Mode of execution                  //
   std::string mode; //Choices are basic [default] and grid [grid manages a fake planet grid]
   
   //Executes a grid of fake planets.
   realT gridCenterSep; ///< The separation of the grid center [pixels].
   realT gridCenterPA; ///< The PA of the grid center [deg E of N].
   realT gridHalfWidthSep; ///< The grid half-width in radius [pixels]
   realT gridDeltaSep; ///< The grid spacing in radius [pixels]
   realT gridHalfWidthPA; ///< The grid half-wdith in PA [pixels]
   realT gridDeltaPA; ///< The grid spacing in PA [pixels]
   std::vector<realT> gridContrasts; ///< The grid contrasts, possibly negative.
   
   int doGrid();
   
   //mx::improc::KLIPreduction<realT, mx::improc::derotVisAO<realT>, realT> * obs;
   mx::improc::KLIPreduction<realT, mx::improc::ADIDerotator<realT>, evCalcT> * obs;
   
   
   
public:
   klipReduce()
   {
      obs = 0;

      deleteFront = 0;
      deleteBack = 0;
      qualityThreshold = 0;
      thresholdOnly = false;
   
      MJDisISO8601 = true;
      imSize = 0;

      angleScale = 0;
      angleConstant = 0;
   
      fakeMethod = 0;
     
      coaddMaxImno = 0;
      coaddMaxTime = 0;
   
      preProcess_beforeCoadd = false;
      preProcess_mask = false;
      preProcess_subradprof = 0;
      preProcess_azUSM_azW = 0;
      preProcess_azUSM_radW = 0;
      preProcess_gaussUSM_fwhm = 0;
      preProcess_only = false;
      skipPreProcess = false;
   
      minDPx = 0;
      includeRefNum = 0;

   
      noDerotate = false;
   
      sigmaThreshold = 0;
   
      exactFName = false;
   
      outputPSFSub = false;
      
      
      
      //Grid setup
      gridCenterSep = 0; 
      gridCenterPA = 0; 
      gridHalfWidthSep = 0; 
      gridDeltaSep = 0; 
      gridHalfWidthPA = 0; 
      gridDeltaPA = 0; 
      
      
      
      
      
      mode = "basic";
   }

   ~klipReduce()
   {
      if(obs) delete obs;
   }
   
   //This sets up the configuration
   void setupConfig()
   {
//       config.add( "help", "h", "help", mx::argType::True,  "", "");
//       config.add( "config","c", "config",mx::argType::Required, "", "config");
//       
      config.add("directory","D", "directory",mx::argType::Required, "", "directory", false, "string", "Directory to search for files");
      config.add("prefix","P", "prefix",mx::argType::Required, "", "prefix", false, "string", "Prefix of the files");
      config.add("extension","E", "extension",mx::argType::Required, "", "extension", false, "string", "");
      config.add("fileList","F", "fileList",mx::argType::Required, "", "fileList", false, "string", "");
      
      config.add("deleteFront","", "deleteFront",mx::argType::Required, "", "deleteFront", false, "int", "");
      config.add("deleteBack","", "deleteBack",mx::argType::Required, "", "deleteBack", false, "int", "");

      config.add("qualityFile","", "qualityFile",mx::argType::Required, "", "qualityFile", false, "string", "");
      config.add("qualityThreshold","", "qualityThreshold",mx::argType::Required, "", "qualityThreshold", false, "", "");
      config.add("thresholdOnly","", "thresholdOnly",mx::argType::True, "", "thresholdOnly", false, "bool", "");
      
      config.add("angleKeyword","", "angleKeyword",mx::argType::Required, "", "angleKeyword", false, "string", "");
      config.add("angleScale","", "angleScale",mx::argType::Required, "", "angleScale", false, "", "");
      config.add("angleConstant","", "angleConstant",mx::argType::Required, "", "angleConstant", false, "", "");
      
      config.add("MJDKeyword","", "MJDKeyword",mx::argType::Required, "", "MJDKeyword", false, "string", "");
      config.add("MJDisISO8601","", "MJDisISO8601",mx::argType::True, "", "MJDisISO8601", false, "", "");
      config.add("imSize","S", "imSize",mx::argType::Required, "", "imSize", false, "", "");
   
      config.add("fakeMethod", "",  "fakeMethod",   mx::argType::Required, "", "fakeMethod", false, "string", "");
      
      config.add("fakeFileName","", "fakeFileName",mx::argType::Required, "", "fakeFileName", false, "string", "");
      config.add("fakeScaleFileName","", "fakeScaleFileName",mx::argType::Required, "", "fakeScaleFileName", false, "string", "");
      config.add("fakeSep","", "fakeSep",mx::argType::Required, "", "fakeSep", false, "", "");
      config.add("fakePA","", "fakePA",mx::argType::Required, "", "fakePA", false, "", "");
      config.add("fakeContrast","", "fakeContrast",mx::argType::Required, "", "fakeContrast", false, "", "");
      
      config.add("coaddMethod","", "coaddMethod",mx::argType::Required, "", "coaddMethod", false, "string", "");
      config.add("coaddMaxImno","", "coaddMaxImno",mx::argType::Required, "", "coaddMaxImno", false, "", "");
      config.add("coaddMaxTime","", "coaddMaxTime",mx::argType::Required, "", "coaddMaxTime", false, "", "");

      config.add("maskFile","", "maskFile",mx::argType::Required, "", "maskFile", false, "string", "");
      
      config.add("preProcess_beforeCoadd","", "preProcess_beforeCoadd",mx::argType::True, "", "preProcess_beforeCoadd", false, "", "");
      config.add("preProcess_mask","", "preProcess_mask",mx::argType::True, "", "preProcess_mask", false, "string", "");
      config.add("preProcess_subradprof","", "preProcess_subradprof",mx::argType::True, "", "preProcess_subradprof", false, "", "");
      config.add("preProcess_azUSM_azW","", "preProcess_azUSM_azW",mx::argType::Required, "", "preProcess_azUSM_azW", false, "", "");
      config.add("preProcess_azUSM_radW","", "preProcess_azUSM_radW",mx::argType::Required, "", "preProcess_azUSM_radW", false, "", "");
      config.add("preProcess_gaussUSM_fwhm","", "preProcess_gaussUSM_fwhm",mx::argType::Required, "", "preProcess_gaussUSM_fwhm", false, "", "");
      config.add("preProcess_outputPrefix","", "preProcess_outputPrefix",mx::argType::Required, "", "preProcess_outputPrefix", false, "", "");
      config.add("preProcess_only","", "preProcess_only",mx::argType::True, "", "preProcess_only", false, "", "");
      config.add("skipPreProcess","", "skipPreProcess",mx::argType::True, "", "skipPreProcess", false, "", "");
      
      config.add("minDPx","", "minDPx",mx::argType::Required, "", "minDPx");
      config.add("excludeMethod","", "excludeMethod",mx::argType::Required, "", "excludeMethod", false, "string", "");
      config.add("includeRefNum","", "includeRefNum",mx::argType::Required, "", "includeRefNum", false, "", "");
      config.add("Nmodes","", "Nmodes",mx::argType::Required, "", "Nmodes", false, "", "");
      config.add("minRadius","", "minRadius",mx::argType::Required, "", "minRadius", false, "", "");
      config.add("maxRadius","", "maxRadius",mx::argType::Required, "", "maxRadius", false, "", "");
      
      config.add("noDerotate","", "noDerotate",mx::argType::False, "", "noDerotate", false, "bool", "Do not derotate before combining.");
      
      config.add("combineMethod",  "", "combineMethod",  mx::argType::Required, "", "combineMethod",  false, "string", "Averaging method for final combination: mean, median, weighted, sigma");
      config.add("weightFile",     "", "weightFile",     mx::argType::Required, "", "weightFile",     false, "string", "File containing weights for the weighted combo.  Two column format: filename weight");
      config.add("sigmaThreshold", "", "sigmaThreshold", mx::argType::Required, "", "sigmaThreshold", false, "float" , "Clipping threshold for sigma clipped mean combination.");
      config.add("minGoodFract", "", "minGoodFract", mx::argType::Required, "", "minGoodFract", false, "float" , "Minimum fraction of good/un-masked pixels to include in final image, otherwise pixel is NaN-ed.");
      
      config.add("outputFile",     "", "outputFile",     mx::argType::Required, "", "outputFile",     false, "string", "Prefix for output file name.  A 4 digit 0-padded number is appended.");
      config.add("exactFName",     "", "exactFName",     mx::argType::True,     "", "exactFName",     false, "bool"  , "Used outputFile exactly as specified, without appending a number or .fits");
      config.add("outputDir",      "", "outputDir",      mx::argType::Required, "", "outputDir",      false, "string", "The directory where to output files.");
      
      config.add("outputPSFSub",      "", "outputPSFSub",      mx::argType::True,     ""    , "outputPSFSub", false, "bool"  , "Output the PSF subtracted images (default false)");
      config.add("psfSubPrefix",      "", "psfSubPrefix",      mx::argType::Required, ""    , "psfSubPrefix", false, "string", "Prefix of the PSF subtracted output files.");
      
      config.add("mode",              "", "mode",              mx::argType::Required, ""    , "mode",         false, "string", "The mode of operation: either \"grid\" or \"normal\" (the default)");
      config.add("grid.centerSep",    "", "grid.centerSep",    mx::argType::Required, "grid", "centerSep",    false, "float" , "The grid center in separation [pixels]" );
      config.add("grid.centerPA",     "", "grid.centerPA",     mx::argType::Required, "grid", "centerPA",     false, "float" , "The grid center in position angle [degrees]" );
      config.add("grid.halfWidthSep", "", "grid.halfWidthSep", mx::argType::Required, "grid", "halfWidthSep", false, "float" , "The half width of the grid in spearation [pixels]" );
      config.add("grid.dalfWidthPA",  "", "grid.halfWidthPA",  mx::argType::Required, "grid", "halfWidthPA",  false, "float" , "The half width of the grid in PA [degrees]" );
      config.add("grid.deltaSep",     "", "grid.deltaSep",     mx::argType::Required, "grid", "deltaSep",     false, "float" , "The grid step size in separation [pixels]" );
      config.add("grid.celtaPA",      "", "grid.deltaPA",      mx::argType::Required, "grid", "deltaPA",      false, "float" , "The grid step size in PA [degrees]" );
      config.add("grid.contrasts",    "", "grid.contrasts",    mx::argType::Required, "grid", "contrasts",    false, "vector float" , "The contrast grid [planet:star]." );
      //config.add("","", "",mx::argType::Required, "", ""));
      
   }

   virtual void setConfigPathCL()
   {
      config.get<std::string>(configPathCL, "config");
   }
   
   void loadConfig()
   { 
      config(directory, "directory");
      config(prefix, "prefix");
      config(extension, "extension");
      config(fileList, "fileList");
      
      config(deleteFront, "deleteFront");
      config(deleteBack, "deleteBack");
      config(qualityFile, "qualityFile");
      config(qualityThreshold, "qualityThreshold");
      config(thresholdOnly, "thresholdOnly");
            
      config(angleKeyword, "angleKeyword");
      config(angleScale,"angleScale");
      config(angleConstant, "angleConstant");
      
      config(MJDKeyword, "MJDKeyword");
      config(MJDisISO8601, "MJDisISO8601");
      config(imSize, "imSize");
      
      config(fakeMethod, "fakeMethod");
      config(fakeFileName, "fakeFileName");
      config(fakeScaleFileName, "fakeScaleFileName");
      
      config(fakeSep, "fakeSep");
      config(fakePA, "fakePA");
      config(fakeContrast, "fakeContrast");
      
      config(coaddMethod, "coaddMethod");
      config(coaddMaxImno, "coaddMaxImno");
      config(coaddMaxTime, "coaddMaxTime");
   
      config(maskFile, "maskFile");
      
      config(preProcess_beforeCoadd, "preProcess_beforeCoadd");
      config(preProcess_mask, "preProcess_mask");
      
      config(preProcess_subradprof, "preProcess_subradprof");
      config(preProcess_azUSM_azW, "preProcess_azUSM_azW");
      config(preProcess_azUSM_radW, "preProcess_azUSM_radW");
      config(preProcess_gaussUSM_fwhm, "preProcess_gaussUSM_fwhm");
      config(preProcess_outputPrefix, "preProcess_outputPrefix");
      config(preProcess_only, "preProcess_only");
      config(skipPreProcess, "skipPreProcess");
      
      config(minDPx, "minDPx");
      config(excludeMethod, "excludeMethod");
      config(includeRefNum, "includeRefNum");
      
      config(Nmodes, "Nmodes");
      
      config(minRadius, "minRadius");
      config(maxRadius, "maxRadius");
      
      config(noDerotate, "noDerotate");
      
      config(combineMethod, "combineMethod");
      
      config(weightFile, "weightFile");
      
      config(sigmaThreshold, "sigmaThreshold");
      config(minGoodFract, "minGoodFract");
      
      config(outputFile, "outputFile");
      config(exactFName, "exactFName");
      
      config(outputDir, "outputDir");
      
      config(outputPSFSub, "outputPSFSub");
      config(psfSubPrefix, "psfSubPrefix");
      
      
      config(gridCenterSep,"grid.centerSep"); 
      
      config(gridCenterPA,"grid.centerPA"); 
      config(gridHalfWidthSep,"grid.halfWidthSep"); 
      config(gridDeltaSep,"grid.deltaSep"); 
      config(gridHalfWidthPA,"grid.halfWidthPA"); 
      config(gridDeltaPA,"grid.deltaPA"); 
      config(gridContrasts, "grid.contrasts");
      
      config(mode, "mode");
      
   }
   
   void printUsage()
   {
      fprintf(stderr, "%s: Perform a KLIP reduction using the mxlib pipeline.\n\n", invokedName.c_str());
      fprintf(stderr, "   usage: %s -D directory -P prefix [-E extension] -n \"x,y,z\" -r X -R x \n\n", invokedName.c_str());
      fprintf(stderr, "   For usage and full documentation see http://makana.as.arizona.edu/acic/group__klipreduce.html  \n\n");
   }
   
   int checkConfig()
   {
      int rv = 0;
      
      if(doHelp) return -1;
      
      
      if(directory == "" && prefix == "" && fileList == "")
      {
         std::cerr << invokedName << ": must specify either directory+prefix or fileList\n\n";
         return -1;
      }
      else
      {
         
         if(fileList != "")
         {
            obs = new mx::improc::KLIPreduction<realT, mx::improc::ADIDerotator<realT>, evCalcT>(fileList);
         }
         else
         {
            obs = new mx::improc::KLIPreduction<realT, mx::improc::ADIDerotator<realT>, evCalcT>(directory, prefix, extension);
         } 
      }
      
      obs->deleteFront = deleteFront;
      obs->deleteBack = deleteBack;
      if(qualityFile != "") obs->qualityFile = qualityFile;
      obs->qualityThreshold = qualityThreshold;
      obs->thresholdOnly = thresholdOnly;
      
      obs->derotF.angleKeyword(angleKeyword);
      obs->derotF.angleScale = angleScale;
      obs->derotF.angleConstant = angleConstant;
      
      if(MJDKeyword != "") obs->MJDKeyword = MJDKeyword;
      obs->MJDisISO8601 = MJDisISO8601;
      obs->imSize = imSize;
      
      obs->fakeMethod = fakeMethod;

      obs->fakeFileName = fakeFileName;
      obs->fakeScaleFileName = fakeScaleFileName;
      obs->fakeSep = fakeSep;
      obs->fakePA = fakePA;
      obs->fakeContrast = fakeContrast;
      
      if(coaddMethod != "")
      {
         if(coaddMethod == "none")
         {
            obs->coaddCombineMethod = mx::improc::HCI::noCombine;
         }
         else if(coaddMethod == "median")
         {
            obs->coaddCombineMethod = mx::improc::HCI::medianCombine;
         }
         else if(coaddMethod == "mean")
         {
            obs->coaddCombineMethod = mx::improc::HCI::meanCombine;
         }
         else
         {
            std::cerr << invokedName << ": invalid coadd method.\n";
            rv = -1;
         }
      }

      obs->coaddMaxImno = coaddMaxImno;
      obs->coaddMaxTime = coaddMaxTime;

      obs->maskFile = maskFile;
      
      obs->preProcess_beforeCoadd = preProcess_beforeCoadd;
      obs->preProcess_mask = preProcess_mask;
      obs->preProcess_subradprof = preProcess_subradprof;
      obs->preProcess_azUSM_azW = preProcess_azUSM_azW;
      obs->preProcess_azUSM_radW = preProcess_azUSM_radW;
      obs->preProcess_gaussUSM_fwhm = preProcess_gaussUSM_fwhm;
      obs->preProcess_outputPrefix = preProcess_outputPrefix;      
      obs->preProcess_only = preProcess_only;
      obs->skipPreProcess = skipPreProcess;
      
      obs->mindpx = minDPx;
      
      if(excludeMethod != "")
      {
         if(excludeMethod == "none")
         {
            obs->excludeMethod = mx::improc::HCI::excludeNone;
         }
         else if(excludeMethod == "pixel")
         {
            obs->excludeMethod = mx::improc::HCI::excludePixel;
         }
         else if(excludeMethod == "angle")
         {
            obs->excludeMethod = mx::improc::HCI::excludeAngle;
         }
         else if(excludeMethod == "imno")
         {
            obs->excludeMethod = mx::improc::HCI::excludeImno;
         }
         else
         {
            std::cerr << invokedName << ": invalid excludeMethod.\n";
            rv = -1;
         }
      }
            
      
      
      
      
      obs->includeRefNum = includeRefNum;
      if(Nmodes.size() == 0)
      {
         std::cerr << invokedName << ": must specify number of modes (Nmodes)\n";
         rv = -1;
      }
      else obs->Nmodes = Nmodes;
      
      if(minRadius.size() == 0 && mode != "postprocess")
      {
         std::cerr << invokedName << ": must specify minimum radii of KLIP regions (minRadius)\n";
         rv = -1;
      }
      
      if(maxRadius.size() == 0 && mode != "postprocess")
      {
         std::cerr << invokedName << ": must specify maximum radii of KLIP regions (maxRadius)\n";
         rv = -1;
      }
      
      obs->doDerotate = !noDerotate;
      
      if(combineMethod != "")
      {
         if(combineMethod == "none")
         {
            obs->combineMethod = mx::improc::HCI::noCombine;
         }
         else if(combineMethod == "median")
         {
            obs->combineMethod = mx::improc::HCI::medianCombine;
         }
         else if(combineMethod == "mean")
         {
            obs->combineMethod = mx::improc::HCI::meanCombine;
         }
         else if(combineMethod == "sigma")
         {
            obs->combineMethod = mx::improc::HCI::sigmaMeanCombine;
         }
         else
         {
            std::cerr << invokedName << ": invalid combine method.\n";
            rv = -1;
         }
      }
      
      obs->weightFile = weightFile;
      obs->sigmaThreshold = sigmaThreshold;
      obs->m_minGoodFract = minGoodFract;
      
      if(outputFile != "")
      {
         obs->finimName = outputFile;
      }
      
      obs->exactFinimName = exactFName;
      
      obs->outputDir = outputDir;
     
      
      obs->doOutputPSFSub = outputPSFSub;
      if(psfSubPrefix != "")
      {
         obs->PSFSubPrefix = psfSubPrefix;
      }
      
      return rv;
   }
   
   

   
   
   virtual int execute()
   {
      if(checkConfig() != 0)
      {
         printUsage();
         return -1;
      }
      
      
      if(mode == "grid")
      {
         return doGrid();
      }
      else if(mode == "postprocess")
      {
         return obs->processPSFSub(directory, prefix, extension);
      }
      else
      {
         std::vector<realT> minMaxQ(minRadius.size(), 0);
         return obs->regions(minRadius, maxRadius, minMaxQ, minMaxQ);
      }
      
   }
   
};

template<typename realT, typename evCalcT>
int klipReduce<realT, evCalcT>::doGrid()
{
   if( gridCenterSep == 0) 
   {
      mxError("klipReduce", MXE_PARAMNOTSET, "Grid center separation not set (grid.centerSep)");
      return -1;
   }
   
   if( gridCenterPA == 0) 
   {
      mxError("klipReduce", MXE_PARAMNOTSET, "Grid center PA not set (grid.centerPA)");
      return -1;
   }
   
   if( gridHalfWidthSep == 0) 
   {
      mxError("klipReduce", MXE_PARAMNOTSET, "Grid half-width in radius not set (grid.halfWidthSep)");
      return -1;
   }
   
   if( gridHalfWidthPA == 0) 
   {
      mxError("klipReduce", MXE_PARAMNOTSET, "Grid half-width in PA not set (grid.halfWidthPA)");
      return -1;
   }
   
   if( gridDeltaSep == 0) 
   {
      mxError("klipReduce", MXE_PARAMNOTSET, "Grid spacing in radius not set (gridDeltaSep)");
      return -1;
   }
   
   if( gridDeltaPA == 0) 
   {
      mxError("klipReduce", MXE_PARAMNOTSET, "Grid spacing in PA not set (grid.deltaPA)");
      return -1;
   }
   
   if( gridContrasts.size() == 0) 
   {
      mxError("klipReduce", MXE_PARAMNOTSET, "Grid contrasts not set (grid.contrasts)");
      return -1;
   }
   
   
   realT x0, y0;
   
   x0 =  -1 * gridCenterSep * sin( mx::math::dtor(gridCenterPA) );
   y0 = gridCenterSep * cos( mx::math::dtor(gridCenterPA) );
   
//   std::cerr << gridCenterSep << " " << x0 << " " << y0 << "\n";
   
   int Nrad = 2 * floor(gridHalfWidthSep / gridDeltaSep) + 1;
   int Npa = 2 * floor(gridHalfWidthPA / gridDeltaPA) + 1;
   
   Eigen::Array<realT, -1, -1> sep, pa;
   
   sep.resize(Nrad, Npa);
   pa.resize(Nrad, Npa);
   
   realT xp, yp,q, x, y;
   
   std::vector<realT> xv, yv;
   
   for(int i=0; i<Nrad; ++i)
   {
      xp = (-0.5*(Nrad - 1) + i)*gridDeltaSep;
      
      for(int j=0; j< Npa; ++j)
      {
         yp = (-0.5*(Npa - 1) + j)*gridDeltaPA;
         
      
         q = mx::math::dtor(90-gridCenterPA);
         
         x = (x0 + xp*cos(q) + yp*sin(q));
         y = (y0 - xp*sin(q) + yp*cos(q));
         
         xv.push_back(x);
         yv.push_back(y);
         
         sep(i,j) = sqrt( pow(x,2) + pow(y,2) );
         
         pa(i,j) = mx::math::angleMod(mx::math::rtod( atan2(y, x))  - 90.0);
         
         //std::cerr << sep(i,j) << " " << pa(i,j) << "\n";
         
         for(size_t k =0; k< gridContrasts.size(); ++k)
         {
            obs->filesRead = false;
            
            obs->fakeSep = {sep(i,j)};
            obs->fakePA = {pa(i,j)};
            obs->fakeContrast = {gridContrasts[k]};
            
            //std::cerr << sep(i,j) << " " << pa(i,j) << " " << gridContrasts[k] << "\n";
            std::vector<realT> minMaxQ(minRadius.size(), 0);
            obs->regions(minRadius, maxRadius, minMaxQ, minMaxQ);
         
         }
         
         
      }
   }
   
   mx::improc::fitsFile<realT> ff;
   
   std::string fn;
   fn = "gridSep.fits";
   if(obs->outputDir != "") fn = outputDir + "/" + fn;
   
   ff.write(fn, sep);
   
   
   fn = "gridPA.fits";
   if(obs->outputDir != "") fn = outputDir + "/" + fn;
   ff.write(fn, pa);
   
   fn = "gridContrasts.dat";
   if(obs->outputDir != "") fn = outputDir + "/" + fn;
   std::ofstream fout;
   fout.open(fn);
   for(size_t i=0; i< gridContrasts.size(); ++i) fout << gridContrasts[i] << "\n";
   fout.close();
   
   return 0;
}
   

int main(int argc, char ** argv)
{
 
#ifdef _OPENMP
   std::cerr << "I think maximum number of threads is: " << omp_get_max_threads() << "\n";
#endif 
   
   
   klipReduce<float, float> kr;
   
   try
   {
      kr.main(argc, argv);
   }
   catch(const mxException & e)
   {
      std::cerr << e.what() << "\n";
      return -1;
   }
   
   return 0;
}

