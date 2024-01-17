/*!
* @file potentialtether.hpp
* @brief ComputeVertexTetherEnergy class
* @author Daniel Matoz, fdamatoz@gmail.com
* @date 20-Sept-2017
*/

#ifndef __potentialtether_HPP__
#define __potentialtether_HPP__

/** @addtogroup computeenergy Compute tethering potential Energy
 *  @brief ComputeVertextetherEnergy definitions
 *  @{
 */
#include <math.h>
//host
#include "computeforceclass.hpp"

//

#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"

/**
 * @class ComputeVertextetherEnergy
 * @brief ComputeVertextetherEnergy just adds a force to vertices of each type
 */
 
/**
 * @brief Example class that acts on the edges 
 * 
 */
class ComputeVertexTetherEnergy : public ComputeForceClass
{
public:
  ComputeVertexTetherEnergy(SystemClass &system) : ComputeForceClass(system)
  {
    m_name = "Tether"; //!< potential name
    m_type = "vertex";  //!< potential type
    this->set_default_properties();
  }
  ~ComputeVertexTetherEnergy() {}

  void set_default_properties(void) override
  {
    pymemb::vector<real> _F(NUM_TYPES_ALLOWED, 0.0);
    m_Fx = _F;
    m_Fy = _F;
    m_Fz = _F;
  }

  using ComputeForceClass::set_property;
  void set_property(std::map<std::string, std::map<std::string, std::string>> &region_map) override
  {
    for (const auto &item : region_map)
    {
      if (item.first.compare("Fx") == 0)
      {
        pymemb::vector<real> _F = util::from_dict_to_vector_types(pymemb::copy(m_Fx), item.second);
        m_Fx = _F;
      }
      else if (item.first.compare("Fy") == 0)
      {
        pymemb::vector<real> _F = util::from_dict_to_vector_types(pymemb::copy(m_Fy), item.second);
        m_Fy = _F;
      }
      else if (item.first.compare("Fz") == 0)
      {
        pymemb::vector<real> _F = util::from_dict_to_vector_types(pymemb::copy(m_Fz), item.second);
        m_Fz = _F;
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
    value["Fx"] = util::to_string_vec(m_Fx);
    value["Fy"] = util::to_string_vec(m_Fy);
    value["Fz"] = util::to_string_vec(m_Fz);
    return value;
  }

  void compute_energy(void) override {}
  void compute(void) override;
  double compute_edge_energy(int) override {return 0.0;}
  double compute_vertex_energy(int) override {return 0.0;}

private:
  pymemb::vector<real> m_Fx; //!< maximum edge allowed
  pymemb::vector<real> m_Fy; //!< maximum edge allowed
  pymemb::vector<real> m_Fz; //!< maximum edge allowed
};

#endif

