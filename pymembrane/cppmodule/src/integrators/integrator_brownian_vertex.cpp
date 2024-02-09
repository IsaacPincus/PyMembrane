#include "integrator_brownian_vertex.hpp"

void IntegratorBrownianMeshVertex::poststep(void)
{
  for (int id = 0; id < _system.Numvertices; id++)
  {
    real3 force_rnd;
    force_rnd.x = force_rnd.y = force_rnd.z = 0.0;
    if (this->get_temperature() > 0.0)
    {
      force_rnd.x = m_B * m_rng->gauss_rng();
      force_rnd.y = m_B * m_rng->gauss_rng();
      force_rnd.z = m_B * m_rng->gauss_rng();
    }
    // Update particle position
    _system.vertices[id].r.x += m_mu * this->get_time_step() * _system.vertices[id].forceC.x + m_sqrt_dt * force_rnd.x;
    _system.vertices[id].r.y += m_mu * this->get_time_step() * _system.vertices[id].forceC.y + m_sqrt_dt * force_rnd.y;
    _system.vertices[id].r.z += m_mu * this->get_time_step() * _system.vertices[id].forceC.z + m_sqrt_dt * force_rnd.z;

    _system.vertices[id].v.x = m_sqrt_dt * force_rnd.x;
    _system.vertices[id].v.y = m_sqrt_dt * force_rnd.y;
    _system.vertices[id].v.z = m_sqrt_dt * force_rnd.z;

    // particles may have been moved slightly outside the box by the above steps, wrap them back into place
    pymemb::enforce_periodic(_system.vertices[id].r, _system.vertices[id].ip, _system.get_box());
  }

  // if there are nanoparticles (not particles as in vertices!!), update them too
  for (int i = 0; i<_system.particles.size(); i++)
  {
    _system.particles[i].position.x += _system.particles[i]._friction*m_mu* this->get_time_step() * _system.particles[i]._total_force.x;
    _system.particles[i].position.y += _system.particles[i]._friction*m_mu* this->get_time_step() * _system.particles[i]._total_force.y;
    _system.particles[i].position.z += _system.particles[i]._friction*m_mu* this->get_time_step() * _system.particles[i]._total_force.z;
  }

}
