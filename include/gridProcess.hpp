

#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <fstream>
#include <cmath>

#include <mx/math/vectorUtils.hpp>
using namespace mx::math;

template<typename realT>
struct gridPt
{
   std::vector<std::string> fnames;
   std::string sep;
   std::string PA;
   std::vector<realT> pas;
   std::string contrast;
   std::string qthresh;
   std::string minrad;
   std::string maxrad;
   std::string mindpx;
   std::string inclrefnum;
   std::string nmodes;
   
   std::vector<realT> snr;
   std::vector<realT> flux;
   std::vector<realT> dr;
   std::vector<realT> dq;
   
   int n;
   
   void setVals(std::string & linein)
   {
      fnames.clear();
      pas.clear();
      snr.clear();
      flux.clear();
      dr.clear();
      dq.clear();
      
      std::istringstream line(linein);
   
      std::string tfn;
      
      line >> tfn;
      fnames.push_back(tfn);
      line >> sep;
      line >> PA;
      
      pas.push_back( mx::ioutils::convertFromString<realT>(PA) );
      
      line >> contrast;
      line >> qthresh;
      line >> minrad;
      line >> maxrad;
      line >> mindpx;
      line >> inclrefnum;
      line >> nmodes;
      
      realT tmp;
      
      line >> tmp;
      snr.push_back(tmp);
      line >> tmp;
      flux.push_back(tmp);
      line >> tmp;
      dr.push_back(tmp);
      line >> tmp;
      dq.push_back(tmp);
      
      
   }
   
   std::string getKey()
   {
      return contrast + "_" + qthresh + "_" + minrad + "_" + maxrad + "_" + mindpx + "_" + inclrefnum +"_" +nmodes;
   }
   
};
   
template<typename realT>
int gridFile( const std::string & outDir,
              const std::string & inName
            )
{
   std::ifstream fin;
   std::string linein;
   gridPt<realT> newpt;
   
   fin.open(inName);

   std::map<std::string, gridPt<realT>> points;
   
   while(!fin.eof())
   {
      getline(fin, linein);
      newpt.setVals(linein);
      newpt.n = 0;
      
      if(points.count(newpt.getKey()) > 0)
      {
         points[newpt.getKey()].fnames.push_back(newpt.fnames[0]);
         points[newpt.getKey()].pas.push_back(newpt.pas[0]);
         points[newpt.getKey()].snr.push_back(newpt.snr[0]);
         points[newpt.getKey()].flux.push_back(newpt.flux[0]);
         points[newpt.getKey()].dr.push_back(newpt.dr[0]);
         points[newpt.getKey()].dq.push_back(newpt.dq[0]);
         points[newpt.getKey()].n++;;
      }
      else
      {
         newpt.n = 1;
         points.insert(std::pair<std::string,gridPt<realT>>(newpt.getKey(), newpt));
      }
   }
   
   
   
   
   typename std::map<std::string, gridPt<realT>>::iterator it = points.begin();

   std::string fn;
   std::ofstream fout;
   while(it != points.end())
   {
      
      fn = outDir + "/" + it->first + ".dat";
      fout.open(fn);
      for(int i=0; i<it->second.n; ++i)
      {
         fout << it->second.fnames[i] << " " << it->second.pas[i] << " " << it->second.flux[i] << " " << it->second.snr[i] << " " << it->second.dr[i] << " " << it->second.dq[i] << "\n";
      }
      fout.close();
      
      
//       if(it->second.contrast != "0.00061") 
//       {
//          ++it;
//          continue;
//       }
      int n = it->second.n;
      std::cout << it->first << " " << it->second.n << " " << it->second.qthresh << " ";
      std::cout << vectorMean(it->second.flux) << " " << sqrt(vectorVariance(it->second.flux)) << " ";
      std::cout << vectorMean(it->second.snr) << " " << sqrt(vectorVariance(it->second.snr)) << "\n";
/*      
      it->second.snr/n << " " << sqrt( it->second.snr_sqr/n - pow(it->second.snr/n,2)) << " ";
      std::cout <<  it->second.flux/n << " " << sqrt( it->second.flux_sqr/n - pow(it->second.flux/n,2)) << " ";
      std::cout << sqrt( it->second.dr_sqr/n - pow(it->second.dr/n,2)) << " ";
      std::cout << sqrt( it->second.dq_sqr/n - pow(it->second.dq/n,2)) <<"\n";*/
      ++it;
   }
   
   
}
