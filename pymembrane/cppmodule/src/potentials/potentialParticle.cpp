#include "potentialParticle.hpp"
#include <cmath>

/*! @fn real ComputeVertexParticleEnergy::compute_energy(int vertex)
    @brief Computes the vertex forces and energies between each vertex in the mesh 
    and the surface of each particle in the system.
    Uses the SDK potential (Soddemann et al. 2001 Eur. Phys. J. E), given as:
    \f$U_{\mathrm{LJ} \cos }= \begin{cases}4\left[\left(\frac{1}{r}\right)^{12}-\left(\frac{1}{r}\right)^6+\frac{1}{4}\right]-\phi, & r \leq 2^{1 / 6}, \\ \frac{1}{2} \phi\left[\cos \left(\alpha r^2+\beta\right)-1\right], & 2^{1 / 6} \leq r \leq 1.5, \\ 0, & r \geq 1.5,\end{cases} \f$
    where \f$\sigma\f$ is the potential range with a well at \f$r = 2^{1/6} \sigma \f$ and zero at \f$r = 1.5 \sigma \f$.
    The strength of the potential is \f$\epsilon \f$, with a well depth of \f$\epsilon \phi \f$.
    @param vertex id
    @return the particle energy at vertex
*/

real get_SDK_energy(real r, 
                    real epsilon,
                    real sigma,
                    real phi)
{
    const real two_sixth_root = 1.12246204831; //!< 2^(1/6)
    const real alpha = 3.1730728678;
    const real beta = -0.85622864544;

    // r must be greater than or equal to 0
    assert(r>=0.0);

    real U;
    if (r<=two_sixth_root*sigma)
    {
        U = 4.0*(pow(sigma/r,12)-pow(sigma/r,6) + 1.0);
    }
    else if (r>=two_sixth_root*sigma && r<=1.5*sigma)
    {
        U = 0.5*phi*(cos(alpha*r*r/(sigma*sigma)+beta)-1);
    }
    else
    {
        U = 0.0;
    }
    return U;
}

real3 get_SDK_force(real d,
                    real3 u_vec, 
                    real epsilon,
                    real sigma,
                    real phi)
{
    const real two_sixth_root = 1.12246204831; //!< 2^(1/6)
    const real alpha = 3.1730728678;
    const real beta = -0.85622864544;

    // r must be greater than or equal to 0
    assert(d>=0.0);

    real Fm;
    real3 F;
    if (d<=two_sixth_root*sigma)
    {
        Fm = 4.0*(-12.0/sigma*pow(sigma/d,13)+6.0/sigma*pow(sigma/d,7));
    }
    else if (d>=two_sixth_root*sigma && d<=1.5*sigma)
    {
        Fm = 0.5*phi*alpha*d/(sigma*sigma)*cos(alpha*d*d/(sigma*sigma)+beta);
    }
    else
    {
        Fm = 0.0;
    }
    Xvec1(F, -Fm, u_vec);
    return F;
}

void ComputeVertexParticleEnergy::compute_energy()
{

    real distance_to_surface;
    real3 ud;
    real U;

    for (int i = 0; i<_particles.size(); i++)
    {
        Particle particle = _particles[i];
        for (int vertex_index = 0; vertex_index < _system.Numvertices; vertex_index++)
        {
            int type = _system.vertices[vertex_index].type;
            real3 r = _system.vertices[vertex_index].r;

            distance_to_surface = particle.get_distance_to_particle_surface(r);

            U = get_SDK_energy(distance_to_surface, this->epsilon[type], this->sigma[type], this->phi[type]);

            _system.vertices[vertex_index].energy += U;
        }
    }
}

void ComputeVertexParticleEnergy::compute() 
{
    real d;
    real3 ud;
    real3 F;
    real3 opp_F;

    for (int i = 0; i<_particles.size(); i++)
    {
        Particle particle = _particles[i];
        for (int vertex_index = 0; vertex_index < _system.Numvertices; vertex_index++)
        {
            int type = _system.vertices[vertex_index].type;
            real3 r = _system.vertices[vertex_index].r;

            d = particle.get_distance_to_particle_surface(r);
            ud = particle.get_unit_vector_to_particle_surface(r);

            F = get_SDK_force(d, ud, this->epsilon[type], this->sigma[type], this->phi[type]);

            _system.vertices[vertex_index].forceC.x += F.x;
            _system.vertices[vertex_index].forceC.y += F.y;
            _system.vertices[vertex_index].forceC.z += F.z;
            // add an equal and opposite force to the particle
            Xvec1(opp_F, -1.0, F);
            _system.particles[i].add_force(r, opp_F);
        }
        
    }
}

real ComputeVertexParticleEnergy::compute_edge_energy(int query_edge_index)
{

}

real ComputeVertexParticleEnergy::compute_vertex_energy(int query_vertex_index)
{

}