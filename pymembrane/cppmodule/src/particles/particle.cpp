#include "particle.hpp"

real Particle::get_distance_to_particle_surface(const real3 &vertex_position)
{
    real3 vector_to_center = pymemb::vector_subtract(vertex_position, this->_position);
    real distance_to_center = sqrt(vdot(vector_to_center,vector_to_center));
    return distance_to_center - this->_radius;
}

real3 Particle::get_unit_vector_to_particle_surface(const real3 &vertex_position)
{
    return pymemb::unit_vector(pymemb::vector_subtract(vertex_position, this->_position));
}

real3 Particle::get_vector_to_particle_surface(const real3 &vertex_position)
{
    return pymemb::vector_subtract(vertex_position, this->_position);
}