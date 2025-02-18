/**
 * \file FMParams.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class flashmatch::FMParams
 *
 * @author kazuhiro
 */

/** \addtogroup core_Base

    @{*/
#ifndef __CVFHICL_FMPARAMS_H__
#define __CVFHICL_FMPARAMS_H__
#include <iostream>
#include <string>
#include <map>
#include "Parser.h"
namespace flashmatch {
  /**
     \class FMParams
     \brief A nested configuration parameter set holder for flashmatch framework.
  */
  class FMParams {
    
  public:
    
    /// Default constructor
    FMParams(const std::string name="",
	 const std::string data="");

    /// Default destructor
    virtual ~FMParams(){};

    /// Copy ctor
    FMParams(const FMParams& orig) : _name       ( orig._name       )
			   , _data_value ( orig._data_value )
			   , _data_pset  ( orig._data_pset  )
    {}

    /// name getter
    inline const std::string& name() const { return _name; }
    
    /// operator override
    inline bool operator==(const FMParams& rhs) const
    {
      if(_name != rhs.name()) return false;
      auto const v_keys = this->value_keys();
      if(v_keys.size() != rhs.value_keys().size()) return false;
      for(auto const& key : v_keys) {
        if(!rhs.contains_value(key))
          return false;
        if(this->get<std::string>(key) != rhs.get<std::string>(key))
          return false;
      }
      auto const p_keys = this->pset_keys();
      if(p_keys.size() != rhs.pset_keys().size()) return false;
      for(auto const& key : p_keys) {
        if(!rhs.contains_pset(key))
          return false;
        if(this->get_pset(key) != rhs.get_pset(key))
          return false;
      }
      return true;
    }

    inline bool operator!=(const FMParams& rhs) const
    { return !((*this) == rhs); }

    /// rename method
    inline void rename(std::string name) { _name = name; }

    /// clear method
    inline void clear() 
    { _data_value.clear(); _data_pset.clear(); }

    /// Set data contents
    void add(const std::string& data);

    /// Insert method for a simple param
    void add_value(std::string key, std::string value);

    /// Insert method for a FMParams rep
    void add_pset(const FMParams& p);

    /// Insert method for a FMParams rep
    void add_pset(std::string key,
		  std::string FMParams);
    
    /// Dump into a text format
    std::string dump(size_t indent_size=0) const;

    /// Dump data string
    std::string data_string() const;

    /// Template getter
    template <class T>
    T get(const std::string& key) const{
      auto iter = _data_value.find(key);
      if( iter == _data_value.end() ) {
        std::string msg;
        msg = "Key does not exist: \"" + key + "\"";
        std::cout<<dump()<<std::endl;
        std::cerr<<msg<<std::endl;
        throw std::exception();
      }
      return parser::FromString<T>((*iter).second);
    }

    /// Template getter w/ default value
    template <class T>
    T get(const std::string& key, const T default_value) const{
      auto iter = _data_value.find(key);
      if( iter == _data_value.end() )
        return default_value;
      return parser::FromString<T>((*iter).second);
    }

    /// None-template function to retrieve parameter set (deprecated)
    const FMParams& get_pset(const std::string& key) const;

    /// Returns # of parameters
    size_t size() const;
    /// Returns a vector of all parameter keys
    const std::vector<std::string> keys() const;
    /// Returns a vector of keys for key-value pairs
    const std::vector<std::string> value_keys () const;
    /// Returns a vector of keys for key-FMParams pairs
    const std::vector<std::string> pset_keys  () const;
    /// Check if a specified key exists for key-value pairs
    bool  contains_value (const std::string& key) const;
    /// Check if a specified key exists for key-FMParams pairs
    bool  contains_pset  (const std::string& key) const;



  private:

    enum KeyChar_t {
      kParamDef,
      kBlockStart,
      kBlockEnd,
      kString,
      kNone
    };

    std::pair<FMParams::KeyChar_t,size_t> search(const std::string& txt, const size_t start) const;
    void strip(std::string& str, const std::string& key);
    void rstrip(std::string& str, const std::string& key);
    void trim_space(std::string& txt);
    void no_space(std::string& txt);

    /// The name of this flashmatch::FMParams
    std::string _name;
    /// Key-Value pairs
    std::map<std::string,std::string> _data_value;
    /// Key-FMParams pairs
    std::map<std::string,::flashmatch::FMParams> _data_pset;

  };

  template<> FMParams FMParams::get<flashmatch::FMParams>(const std::string& key) const;
  
}

#endif
/** @} */ // end of doxygen group 

