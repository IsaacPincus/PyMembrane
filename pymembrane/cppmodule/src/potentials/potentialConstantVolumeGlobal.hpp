#ifndef __potentialConstantVolumeGlobal_HPP__
#define __potentialConstantVolumeGlobal_HPP__

#include "computeforceclass.hpp"
#include "../system/systemclass.hpp"
#include "../mesh/computegeometry.hpp"

/**
 * @class ComputeVertexConstantGlobalVolumeEnergy
 * @brief ComputeVertexConstantGlobalVolumeEnergy Compute Vertex Volume Gradient
 *  For a given triangle t compute contribution from the surface tension caused by change in the global surface area
 *
 *     \f$ E = \kappa_{v} \frac{(V-V_0)^2}{V_0} \f$
 *
 *     where \f$ V \f$ is current global volume, \f$ V_0 \f$  is the native (given) global volume, and
 *     \f$ \kappa_{v} \f$ is the volume expansion modulus.
 * 
 *     Note that this follows Bian et al. 2020 in terms of nomenclature.
 *
 */

class ComputeVertexConstantGlobalVolumeEnergy : public ComputeForceClass
{
public:
    ComputeVertexConstantGlobalVolumeEnergy(SystemClass &system) : ComputeForceClass(system)
    {
        m_name = "Constant Global Volume"; //!< potential name
        m_type = "global";          //!< potential type
        this->set_default_properties();
    }

    ~ComputeVertexConstantGlobalVolumeEnergy() {}

    void set_default_properties(void) override
    {
        real _kappa_v;
        kappa_v = _kappa_v;
        real _target_volume;
        target_volume = _target_volume;

        // get global volume of this membrane
        current_volume = 0;
        if (_system.close_surface == true)
        {
            for (int face_index = 0; face_index < _system.faces.size(); face_index++)
            {
                real3 r0 = _system.vertices[_system.faces[face_index].v1].r;
                real3 r1 = _system.vertices[_system.faces[face_index].v2].r;
                real3 r2 = _system.vertices[_system.faces[face_index].v3].r;
                real3 nt = pymemb::compute_normal_triangle(r0, r1, r2, _system.get_box());
                current_volume += vdot(r0, nt) / 6.0;
            }
        }
        else
        {
            std::cerr << "open surface doesn't have any volume\n";
        }
    }

    using ComputeForceClass::set_property;
    void set_property(std::map<std::string, std::string> &region_map) override
    {
        for (const auto &item : region_map)
        {
            auto propname = item.first;
            auto value = item.second;
            if (propname.compare("kappa_v") == 0)
            {
                real _kappa_v = util::from_string_double(value);
                kappa_v = _kappa_v;
            }
            else if (propname.compare("target_volume") == 0)
            {
                real _target_volume = util::from_string_double(value);
                target_volume = _target_volume;
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
        value["kappa_v"] = util::to_string(kappa_v);
        value["target_volume"] = util::to_string(target_volume);

        return value;
    }

    void compute_energy(void) override;
    void compute(void) override;
    real compute_edge_energy(int) override;
    real compute_vertex_energy(int) override;
    //void compute_stress(void) override;
    //void compute_atomic_stress(void) override;

private:
    real kappa_v;       //!< strength
    real target_volume; //!< target global equilibirum area
    real current_volume;
    void recalculate_mesh_volume(void);
};

#endif
