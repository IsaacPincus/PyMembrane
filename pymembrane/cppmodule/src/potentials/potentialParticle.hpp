/*!
* @file potentialParticle.hpp
* @brief ComputeVertexParticleEnergy class
* @author Isaac Pincus, isaac.pincus@gmail.com
* @date 12-Jan-2024
*/
#ifndef __POTENTIALPARTICLE_HPP__
#define __POTENTIALPARTICLE_HPP__

/** @addtogroup computeenergy Compute Vertex Particle Energy
 *  @brief ComputeVertexParticleEnergy definitions
 *  @{
 */
//host
#include "computeforceclass.hpp"

//
#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"
#include "../particles/particle.hpp"

/**
 * @class ComputeVertexParticleEnergy
 * @brief ComputeVertexParticleEnergy Compute the Particle-Membrane Potential
 */
class ComputeVertexParticleEnergy : public ComputeForceClass
{
public:
    ComputeVertexParticleEnergy(SystemClass &system) : ComputeForceClass(system)
    {
        m_name = "Particle"; //!< potential name
        m_type = "vertex";    //!< potential type
        _particles = _system.particles;
        this->set_default_properties();
    }
    
    void set_default_properties(void) override
    {
        pymemb::vector<real> _epsilon(NUM_TYPES_ALLOWED, 0.0);
        epsilon = _epsilon;
        pymemb::vector<real> _sigma(NUM_TYPES_ALLOWED, 1.0);
        sigma = _sigma;
        pymemb::vector<real> _phi(NUM_TYPES_ALLOWED, 1.0);
        phi = _phi;
    }

    using ComputeForceClass::set_property;
    void set_property(std::map<std::string, std::map<std::string, std::string>> &region_map) override
    {
        for (const auto &item : region_map)
        {
        if (item.first.compare("epsilon") == 0)
        {
            pymemb::vector<real> _epsilon = util::from_dict_to_vector_types(pymemb::copy(epsilon), item.second);
            epsilon = _epsilon;
        }
        else if (item.first.compare("sigma") == 0)
        {
            pymemb::vector<real> _sigma = util::from_dict_to_vector_types(pymemb::copy(sigma), item.second);
            sigma = _sigma;
        }
        else if (item.first.compare("phi") == 0)
        {
            pymemb::vector<real> _phi = util::from_dict_to_vector_types(pymemb::copy(phi), item.second);
            phi = _phi;
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
        value["epsilon"] = util::to_string_vec(epsilon);
        value["sigma"] = util::to_string_vec(sigma);
        value["phi"] = util::to_string_vec(phi);
        return value;
    }

  void compute_energy(void) override;
  void compute(void) override;
  real compute_edge_energy(int) override;
  real compute_vertex_energy(int) override;

private:
    pymemb::vector<real> epsilon; //!< bending rigidity
    pymemb::vector<real> sigma; //!< bending rigidity
    pymemb::vector<real> phi; //!< bending rigidity
    pymemb::vector<Particle> _particles; // pass this by reference so we can update it
};

#endif