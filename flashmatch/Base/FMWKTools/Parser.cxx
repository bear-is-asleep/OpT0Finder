#ifndef __FLASHMATCHBASE_PARSER_CXX__
#define __FLASHMATCHBASE_PARSER_CXX__

#include "PSetUtils.h"
namespace flashmatch {
  namespace parser{
    
    template<> std::string FromString( const std::string& value)
    {
      std::string res(value);
      if(res.empty()) return res;
      
      if(res.find("\"") == 0) res = res.substr(1);
      if(res.empty()) return res;
      
      if(res.rfind("\"") == (res.length()-1)) res = res.substr(0,res.length()-1);
      return res;
    }
    
    template<> float FromString( const std::string& value )
    { return std::stof(value); }
    
    template<> double FromString( const std::string& value )
    { return std::stod(value); }
    
    template<> unsigned short FromString( const std::string& value)
    { return std::stoul(value); }
    
    template<> unsigned int FromString( const std::string& value)
    { return std::stoul(value); }
    
    template<> unsigned long FromString( const std::string& value)
    { return std::stoul(value); }

    template<> short FromString( const std::string& value )
    { return std::stoi(value); }
    
    template<> int FromString( const std::string& value )
    { return std::stoi(value); }
    
    template<> long FromString( const std::string& value )
    { return std::stol(value); }
    
    template<> bool FromString(const std::string& value )
    {
      std::string tmp=value;
      std::transform(tmp.begin(),tmp.end(),tmp.begin(),::tolower);
      if( value == "true"  || value == "1" ) return true;
      if( value == "false" || value == "0" ) return false;
      std::cerr << "Invalid bool expression: " << tmp << std::endl;
      throw std::exception();
      return false;
    }
    
    template<> std::vector< std::string    > FromString (const std::string& value )
    {
      std::vector<std::string> res;
      if(value.find("[") != 0 || (value.rfind("]")+1) != value.size()) {
        std::string msg;
	      std::cerr << "Invalid vector expression: " << value << std::endl;
        throw std::exception();
      }
      size_t index = 1;
      while((index+1) < value.size()) {
        size_t next_index = value.find(",",index);
        if(next_index >= value.size()) break;
        std::string cand = value.substr(index,next_index-index);
        res.emplace_back(cand);
        index = next_index + 1;
      }
      if( (index+1) < value.size() )
        res.push_back(value.substr(index,value.size()-index-1));

      for(auto& s : res) {
        if(s.find("\"")==0) s=s.substr(1,s.size()-1);
        if(s.rfind("\"")+1 == s.size()) s = s.substr(0,s.size()-1);
      }
      return res;
    }

    template<> std::vector< std::vector< std::string > > FromString( const std::string& value )
    {
      //Use split2d and trim
      std::vector<std::vector<std::string>> res;
      std::string trimmed = Trim(value);
      std::vector<std::vector<std::string>> split2d = Split2D(trimmed);
      for (size_t i = 0; i < split2d.size(); i++) {
        std::vector<std::string> tmp;
        for (size_t j = 0; j < split2d[i].size(); j++) {
          tmp.push_back(split2d[i][j]);
        }
        res.push_back(tmp);
      }
      return res;
    }
    
    template<> std::vector< float > FromString< std::vector< float > > (const std::string& value )
    {
      auto str_v = FromString<std::vector<std::string> >(value);
      std::vector<float> res;
      res.reserve(str_v.size());
      for(auto const& v : str_v)
        res.push_back( FromString<float>(v) );
      return res;
    }
    
    template<> std::vector< double > FromString< std::vector< double > > (const std::string& value )
    {
      auto str_v = FromString<std::vector<std::string> >(value);
      std::vector<double> res;
      res.reserve(str_v.size());
      for(auto const& v : str_v)
        res.push_back( FromString<double>(v) );
      return res;
    }
    
    template<> std::vector< short > FromString< std::vector< short > > (const std::string& value )
    {
      auto str_v = FromString<std::vector<std::string> >(value);
      std::vector<short> res;
      res.reserve(str_v.size());
      for(auto const& v : str_v)
        res.push_back( FromString<short>(v) );
      return res;
    }
    
    template<> std::vector< int > FromString< std::vector< int > > (const std::string& value )
    {
      auto str_v = FromString<std::vector<std::string> >(value);
      std::vector<int> res;
      res.reserve(str_v.size());
      for(auto const& v : str_v)
        res.push_back( FromString<int>(v) );
      return res;
    }
    
    template<> std::vector< long > FromString< std::vector< long > > (const std::string& value )
    {
      auto str_v = FromString<std::vector<std::string> >(value);
      std::vector<long> res;
      res.reserve(str_v.size());
      for(auto const& v : str_v)
        res.push_back( FromString<long>(v) );
      return res;
    }
    
    template<> std::vector< unsigned short > FromString< std::vector< unsigned short > > (const std::string& value )
    {
      auto str_v = FromString<std::vector<std::string> >(value);
      std::vector<unsigned short> res;
      res.reserve(str_v.size());
      for(auto const& v : str_v)
        res.push_back( FromString<unsigned short>(v) );
      return res;
    }
    
    template<> std::vector< unsigned int > FromString< std::vector< unsigned int > > (const std::string& value )
    {
      auto str_v = FromString<std::vector<std::string> >(value);
      std::vector<unsigned int> res;
      res.reserve(str_v.size());
      for(auto const& v : str_v)
        res.push_back( FromString<unsigned int>(v) );
      return res;
    }
    
    template<> std::vector< unsigned long > FromString< std::vector< unsigned long > > (const std::string& value )
    {
      auto str_v = FromString<std::vector<std::string> >(value);
      std::vector<unsigned long> res;
      res.reserve(str_v.size());
      for(auto const& v : str_v)
        res.push_back( FromString<unsigned long>(v) );
      return res;
    }

    template<> std::vector< bool > FromString< std::vector< bool > > (const std::string& value )
    {
      auto str_v = FromString<std::vector<std::string> >(value);
      std::vector<bool> res;
      res.reserve(str_v.size());
      for(auto const& v : str_v)
        res.push_back( FromString<bool>(v) );
      return res;
    }

    template<> std::vector<std::vector<double>> FromString<std::vector<std::vector<double>>>(const std::string& value) {
        auto str_v = FromString<std::vector<std::string>>(value);
        std::vector<std::vector<double>> res;
        std::string trimmed = Trim(value);
        std::vector<std::string> split = Split(trimmed);
        for(auto const& v : split) {
            std::vector<double> tmp = FromString<std::vector<double>>(v);
            res.push_back(tmp);
        }
        return res;
    }

    template<> std::vector<std::vector<std::vector<double>>> FromString<std::vector<std::vector<std::vector<double>>>>(const std::string& value) {
        auto str_v = FromString<std::vector<std::vector<std::string>>>(value);
        std::vector<std::vector<std::vector<double>>> res(str_v.size());
        for (size_t i = 0; i < str_v.size(); i++) {
            for (size_t j = 0; j < str_v[i].size(); j++) {
              std::vector<double> tmp = FromString<std::vector<double>>(str_v[i][j]);
              res[i].push_back(tmp);
            }
        }
        return res;
    }
    
    template<> std::string ToString<std::string>(const std::string& value)
    {
      std::string res(value);
      if(res.empty()) return res;
      
      if(res.find("\"") == 0) res = res.substr(1);
      if(res.empty()) return res;
      
      if(res.rfind("\"") == (res.length()-1)) res = res.substr(0,res.length()-1);
      return res;
    }

    std::string Trim(const std::string& value)
    {
      std::string res(value);
      if(res.empty()) return res;
      
      if(res.find("\"") == 0) res = res.substr(1);
      if(res.empty()) return res;
      
      if(res.rfind("\"") == (res.length()-1)) res = res.substr(0,res.length()-1);
      while (res.find("\n") != std::string::npos) {
        res.erase(res.find("\n"), 1);
      }
      while (res.find("\t") != std::string::npos) {
        res.erase(res.find("\t"), 1);
      }
      while (res.find(" ") != std::string::npos) {
        res.erase(res.find(" "), 1);
      }

      return res;
    }

    std::vector<std::string> Split(const std::string& value, const std::string& delim)
    {
      std::vector<std::string> res;
      if(value.empty()) return res;
      size_t index = 0;
      while(index < value.size()) {
        size_t next_index = value.find(delim,index);
        if(next_index == std::string::npos) break;
        std::string cand = value.substr(index,next_index-index);
        //Remove all [] from the string
        while (cand.find("[") != std::string::npos) {
          cand.erase(cand.find("["), 1);
        }
        while (cand.find("]") != std::string::npos) {
          cand.erase(cand.find("]"), 1);
        }
        //Break if it's empty
        if (cand.empty()) break;
        if (cand == delim) break;
        if (std::all_of(cand.begin(),cand.end(),::isspace)) break; //all spaces, no values
        //Remove leading commas
        if (cand.find(",") == 0) {
          cand = cand.substr(1, cand.size() - 1);
        }
        //Add a [] to edges to convert to vector
        cand = "[" + cand + "]";
        res.emplace_back(cand);
        index = next_index + delim.size();
      }
      
      return res;
    }

    std::vector<std::vector<std::string>> Split2D(const std::string& value, const std::string& delim)
    {
      std::vector<std::vector<std::string>> res;
      if(value.empty()) return res;
      size_t index = 0;
      while(index < value.size()) {
        size_t next_index = value.find(delim,index);
        if(next_index == std::string::npos) break;
        std::string cand = value.substr(index,next_index-index);
        //Remove all [] from the string
        while (cand.find("[[") != std::string::npos) {
          cand.replace(cand.find("[["), 2, "[");
        }
        while (cand.find("]]") != std::string::npos) {
          cand.replace(cand.find("]]"), 2, "]");
        }
        //Break if it's empty
        if (cand.empty()) break;
        if (cand == delim) break;
        if (std::all_of(cand.begin(),cand.end(),::isspace)) break; //all spaces, no values
        //Remove leading commas
        if (cand.find(",") == 0) {
          cand = cand.substr(1, cand.size() - 1);
        }
        //Add a [[]] to edges to convert to matrix
        cand = "[[" + cand + "]]";
        res.emplace_back(Split(cand));
        index = next_index + delim.size();
      }
      
      return res;
    }
    
  }
}
#endif
  
