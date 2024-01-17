#ifndef __POTENTIALFENEPOW_HPP__
#define __POTENTIALFENEPOW_HPP__


#include "computeforceclass.hpp"
#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"

/**
 *  @class ComputeVertexFENEPOWEnergy
 *  @brief Implements methods to compute the stretching contribution to the edge energy.
 * 
 *  Calculates stretching energy of the edge using:
 *  \f$ E_{FENE-POW} = \frac{k_s}{2} l_\mathrm{max}^2 \log \left[ 1-x^2 \right] \f$ where \f$ k_s \f$ is the spring constant,
 *  \f$ l \f$ is current edge length, l_\mathrm{max} is the maximum edge length, and \f$ x = l/l_\mathrm{max} \f$.
 *  There is also a repulsive power function potential with the form \f$ k_p/l^m \f$.
 */

class ComputeVertexFENEPOWEnergy : public ComputeForceClass
{
public:
  ComputeVertexFENEPOWEnergy(SystemClass &system) : ComputeForceClass(system)
  {
    m_name = "FENEPOW"; //!< potential name
    m_type = "edge";     //!< potential type
    this->set_default_properties();
  }

  ~ComputeVertexFENEPOWEnergy() {}

  void set_default_properties(void) override
  {
    pymemb::vector<real> _k(NUM_TYPES_ALLOWED, 0.0);
    m_ks = _k;
    pymemb::vector<real> _lmax(NUM_TYPES_ALLOWED, 1.0);
    m_lmax = _lmax;
    pymemb::vector<real> _kp(NUM_TYPES_ALLOWED, 1.0);
    m_kp = _kp;
    pymemb::vector<real> _m(NUM_TYPES_ALLOWED, 2.0);
    m_m = _m;
  }

  using ComputeForceClass::set_property;
  void set_property(std::map<std::string, std::map<std::string, std::string>> &region_map) override
  {
    for (const auto &item : region_map)
    {
      if (item.first.compare("ks") == 0)
      {
        pymemb::vector<real> _ks = util::from_dict_to_vector_types(pymemb::copy(m_ks), item.second);
        m_ks = _ks;
      }
      else if (item.first.compare("lmax") == 0)
      {
        pymemb::vector<real> _lmax = util::from_dict_to_vector_types(pymemb::copy(m_lmax), item.second);
        m_lmax = _lmax;
      }
      else if (item.first.compare("kp") == 0)
      {
        pymemb::vector<real> _kp = util::from_dict_to_vector_types(pymemb::copy(m_kp), item.second);
        m_kp = _kp;
      }
      else if (item.first.compare("m") == 0)
      {
        pymemb::vector<real> _ml = util::from_dict_to_vector_types(pymemb::copy(m_m), item.second);
        m_m = _ml;
        // check all m values are greater than 0
        for (real element : m_m)
        {
          assert(element > 1);
        }
      }
      else
        this->print_warning_property_name(item.first);
    }
  }
  
  std::map<std::string, std::string> get_info(void) override
  {
    std::map<std::string, std::string> value;
    value["name"] = m_name;
    value["type"] = m_type;
    value["ks"] = util::to_string_vec(m_ks);
    value["lmax"] = util::to_string_vec(m_lmax);
    value["kp"] = util::to_string_vec(m_kp);
    value["m"] = util::to_string_vec(m_m);
    return value;
  }

  void compute_energy(void) override;
  void compute(void) override;
  real compute_edge_energy(int) override;
  real compute_vertex_energy(int) override;
  void compute_stress(void) override;
  void compute_atomic_stress(void) override;

private:
  pymemb::vector<real> m_ks;  //!< spring constant
  pymemb::vector<real> m_lmax; //!< maximum spring length
  pymemb::vector<real> m_kp; //!< power law spring constant
  pymemb::vector<real> m_m; //!< power law exponent
};

#endif

