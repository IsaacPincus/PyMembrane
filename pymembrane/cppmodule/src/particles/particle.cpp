#include "particle.hpp"

real Particle::get_distance_to_particle_surface(const real3 &vertex_position)
{
    real3 vector_to_center = pymemb::vector_subtract(vertex_position, this->position);
    real distance_to_center = sqrt(vdot(vector_to_center,vector_to_center));
    return distance_to_center - this->_radius;
}

real3 Particle::get_unit_vector_to_particle_surface(const real3 &vertex_position)
{
    return pymemb::unit_vector(pymemb::vector_subtract(vertex_position, this->position));
}

real3 Particle::get_vector_to_particle_surface(const real3 &vertex_position)
{
    real distance_to_surface = this->get_distance_to_particle_surface(vertex_position);
    real3 unit_vec_to_surface = this->get_unit_vector_to_particle_surface(vertex_position);
    real3 vec_to_surface;
    vec_to_surface.x = unit_vec_to_surface.x*distance_to_surface;
    vec_to_surface.y = unit_vec_to_surface.y*distance_to_surface;
    vec_to_surface.z = unit_vec_to_surface.z*distance_to_surface;
    return vec_to_surface;
}

real3 Particle::get_closest_point_on_surface(const real3 &vertex_position)
{
    real3 unit_vec_to_surface = this->get_unit_vector_to_particle_surface(vertex_position);
    // reverse unit vector and then multiply by radius
    real3 point_on_surface;
    point_on_surface.x = -unit_vec_to_surface.x*this->_radius;
    point_on_surface.y = -unit_vec_to_surface.y*this->_radius;
    point_on_surface.z = -unit_vec_to_surface.z*this->_radius;
}

void Particle::add_force(const real3 &vertex_position, const real3 &force)
{
    this->_total_force.x += force.x;
    this->_total_force.y += force.y;
    this->_total_force.z += force.z;

    // in the future, may want to also update the particle torque, requiring the point where the force is applied
    // real3 closest_point_on_surface = this->get_closest_point_on_surface(vertex_position);
}