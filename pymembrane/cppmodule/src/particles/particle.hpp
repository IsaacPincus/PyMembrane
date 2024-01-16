#ifndef __particle_hpp__
#define __particle_hpp__

#pragma once

#include "../types/globaltypes.hpp"
#include "../mesh/computegeometry.hpp"

// for now we only have spherical particles
class Particle {
    public:
        Particle() 
        {
            this->set_default_properties();
        }
        Particle(real3 position, real radius) : _position(position), _radius(radius)
        {
        }
        Particle(real x, real y, real z, real radius) : _radius(radius)
        {
            _position.x = x;
            _position.y = y;
            _position.z = z;
        }

        void set_default_properties(void)
        {
            _radius = 1.0;
            _position.x = 0.0; _position.y = 0.0; _position.z = 0.0;
        }

        real get_distance_to_particle_surface(const real3 &vertex_position);
        real3 get_unit_vector_to_particle_surface(const real3 &vertex_position);
        real3 get_vector_to_particle_surface(const real3 &vertex_position);

    private:
        real _radius;
        real3 _position;
};

#endif