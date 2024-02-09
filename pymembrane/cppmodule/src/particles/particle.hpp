#ifndef __particle_hpp__
#define __particle_hpp__

#pragma once

#include "../types/globaltypes.hpp"
#include "../mesh/computegeometry.hpp"

#include <pybind11/pybind11.h>
namespace py = pybind11;

// for now we only have spherical particles
class Particle {
    public:
        Particle() 
        {
            this->set_default_properties();
        }
        Particle(real3 position, real radius) : position(position), _radius(radius)
        {
        }
        Particle(real x, real y, real z, real radius) : _radius(radius)
        {
            position.x = x;
            position.y = y;
            position.z = z;
            _total_force.x = 0.0;
            _total_force.y = 0.0;
            _total_force.z = 0.0;
            _friction = 1.0/(6.0*defPI*_radius);
        }

        void set_default_properties(void)
        {
            _radius = 1.0;
            position.x = 0.0; position.y = 0.0; position.z = 0.0;
        }

        void reset_forces() 
        {
            _total_force.x = 0.0;
            _total_force.y = 0.0;
            _total_force.z = 0.0;
        }

        real get_distance_to_particle_surface(const real3 &vertex_position);
        real3 get_unit_vector_to_particle_surface(const real3 &vertex_position);
        real3 get_vector_to_particle_surface(const real3 &vertex_position);
        real3 get_closest_point_on_surface(const real3 &vertex_position);
        void add_force(const real3 &vertex_position, const real3 &force);


        real3 position;
        real3 _total_force;
        real _friction;

    private:
        real _radius;
        
};

#endif