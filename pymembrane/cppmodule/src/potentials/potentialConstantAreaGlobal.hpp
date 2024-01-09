#ifndef __potentialConstantAreaGlobal_HPP__
#define __potentialConstantAreaGlobal_HPP__

#include "computeforceclass.hpp"
#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"

/**
 * @class ComputeVertexConstantGlobalAreaEnergy
 * @brief ComputeVertexConstantGlobalAreaEnergy Compute Vertex Cauchy-Green Stretching Energy
 *  For a given triangle t compute contribution from the surface tension caused by change in the global surface area
 *
 *     \f$ E = \kappa_{ag} \frac{(A-A_0)^2}{A_0} \f$
 *
 *     where \f$ A \f$ is current global surface area, \f$ A_0 \f$  is the native (given) global surface area, and
 *     \f$ \kappa_{ag} \f$ is the stretching modulus.
 * 
 *     Note that this follows Bian et al. 2020 in terms of nomenclature.
 *
 */

class ComputeVertexConstantGlobalAreaEnergy : public ComputeForceClass
{
public:
    ComputeVertexConstantGlobalAreaEnergy(SystemClass &system) : ComputeForceClass(system)
    {
        m_name = "Constant Global Area"; //!< potential name
        m_type = "global";          //!< potential type
        this->set_default_properties();
    }

    ~ComputeVertexConstantGlobalAreaEnergy() {}

    void set_default_properties(void) override
    {
        real _kappa_ag;
        kappa_ag = _kappa_ag;
        real _target_area;
        target_area = _target_area;

        // get global area of this membrane
        current_area = 0;
        for (int face_index = 0;
            face_index < _system.Numfaces;
            face_index++)
        {
            int v1 = _system.faces[face_index].v1;
            int v2 = _system.faces[face_index].v2;
            int v3 = _system.faces[face_index].v3;
            current_area += pymemb::compute_area_triangle_from_vertex(
                    _system.vertices[v1].r,_system.vertices[v2].r, _system.vertices[v3].r, _system.get_box());
        }
        
    }

    using ComputeForceClass::set_property;
    void set_property(std::map<std::string, std::string> &region_map) override
    {
        for (const auto &item : region_map)
        {
            auto propname = item.first;
            auto value = item.second;
            if (propname.compare("kappa_ag") == 0)
            {
                real _kappa_ag = util::from_string_double(value);
                kappa_ag = _kappa_ag;
            }
            else if (propname.compare("target_area") == 0)
            {
                real _target_area = util::from_string_double(value);
                target_area = _target_area;
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
        value["kappa_ag"] = util::to_string(kappa_ag);
        value["target_area"] = util::to_string(target_area);

        return value;
    }

    void compute_energy(void) override;
    void compute(void) override;
    real compute_edge_energy(int) override;
    real compute_vertex_energy(int) override;
    // void compute_stress(void) override;
    // void compute_atomic_stress(void) override;

private:
    real kappa_ag;       //!< strength
    real target_area; //!< target global equilibirum area
    real current_area;
    void recalculate_mesh_area(void);
};

#endif
