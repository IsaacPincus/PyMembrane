#include "potentialConstantVolumeGlobal.hpp"

void ComputeVertexConstantGlobalVolumeEnergy::recalculate_mesh_volume(void)
{
    // get global volume of this membrane
    current_volume = 0;
    if (_system.close_surface == true)
    {
        for (int face_index = 0; face_index < _system.faces.size(); face_index++)
        {
            // see Bian et al. 2020 Appendix G. Orientation of surface (signed volume) takes care of interior/exterior.
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

real ComputeVertexConstantGlobalVolumeEnergy_lambda(const real _current_volume,
                                                    const real _kappa_v,
                                                    const real _target_volume)
{
    real energy = _kappa_v/_target_volume * (_current_volume - _target_volume) * (_current_volume - _target_volume);
    return energy;
}

void ComputeVertexConstantGlobalVolumeEnergy_fn(const int Numfaces,
                                                HE_FaceProp *faces,
                                                const real _current_volume,
                                                const real _kappa_v,
                                                const real _target_volume)
{

    real energy = ComputeVertexConstantGlobalVolumeEnergy_lambda(_current_volume, _kappa_v, _target_volume);
    for (int face_index = 0; face_index < Numfaces; face_index++)
    {
        /// Add energy to that face
        faces[face_index].energy += energy;
    }
}

void ComputeVertexConstantGlobalVolumeEnergy::compute_energy(void)
{

    ComputeVertexConstantGlobalVolumeEnergy_fn(_system.Numfaces,
                                               &_system.faces[0],
                                               current_volume,
                                               kappa_v,
                                               target_volume);
}

real ComputeVertexConstantGlobalVolumeEnergy::compute_edge_energy(int query_edge_index)
{
    // reset area to recalculate
    this->recalculate_mesh_volume();
    // we need to loop the 2 faces that are connected to the edge_index
    pymemb::vector<int> face_vec{_system.edges[query_edge_index].face_k, _system.edges[query_edge_index].face_l};
    // reset energy
    real edge_energy = 0.0;
    for (auto face_index : face_vec)
    {
        edge_energy += ComputeVertexConstantGlobalVolumeEnergy_lambda(current_volume, kappa_v, target_volume);
    }
    return edge_energy;
}

real ComputeVertexConstantGlobalVolumeEnergy::compute_vertex_energy(int query_vertex_index)
{
    // reset area to recalculate
    this->recalculate_mesh_volume();
    real energy = 0.0;

    ///< get the triangle that this vertex is part of
    int he = _system.vertices[query_vertex_index]._hedge;
    int first = he;
    // std::cout<< "first " << first << "\n";
    int face_index, he_pair, he_pair_next;
    do
    {

        face_index = _system.halfedges[he].face;
        if (_system.faces[face_index].boundary == false) // Remember -1 is the virtual face outside of the mesh
        {
            energy += kappa_v*(current_volume-target_volume)*(current_volume-target_volume)/(target_volume);
        }

        int he_prev = _system.halfedges[he].prev;
        he = _system.halfedges[he_prev].pair;
    } while ((he != first));
    return energy;
}

void ComputeVertexConstantVolumeTriangleForce_fn(const int Numfaces,
                                               HE_VertexProp *vertices,
                                               const HE_FaceProp *faces,
                                               const real _kappa_v,
                                               const real _target_volume,
                                               const real _current_volume,
                                               const BoxType _box)
{

    real current_volume = 0;
    real force_factor = -_kappa_v*2*(_current_volume-_target_volume)/(_target_volume);
    
    for (int face_index = 0; face_index < Numfaces; face_index++)
    {
        int type = faces[face_index].type;
        int v1 = faces[face_index].v1;
        int v2 = faces[face_index].v2;
        int v3 = faces[face_index].v3;

        real3 r1 = vertices[v1].r;
        real3 r2 = vertices[v2].r;
        real3 r3 = vertices[v3].r;

        /*----------------------------------------------------------------------------------------------------------------*/
        /*-----------------------------------           ACTUAL CALCULATION        ----------------------------------------*/
        /*----------------------------------------------------------------------------------------------------------------*/
        real3 cross;
        // v1
        cross = pymemb::vector_cross(r2, r3);
        vertices[v1].forceC.x += 1/6*force_factor * cross.x;
        vertices[v1].forceC.y += 1/6*force_factor * cross.y;
        vertices[v1].forceC.z += 1/6*force_factor * cross.z;

        // v2
        cross = pymemb::vector_cross(r1, r3);
        vertices[v2].forceC.x += 1/6*force_factor * cross.x;
        vertices[v2].forceC.y += 1/6*force_factor * cross.y;
        vertices[v2].forceC.z += 1/6*force_factor * cross.z;

        // v3
        cross = pymemb::vector_cross(r1, r2);
        vertices[v3].forceC.x += 1/6*force_factor * cross.x;
        vertices[v3].forceC.y += 1/6*force_factor * cross.y;
        vertices[v3].forceC.z += 1/6*force_factor * cross.z;
    }
}

void ComputeVertexConstantGlobalVolumeEnergy::compute(void)
{

    this->recalculate_mesh_volume();

    ComputeVertexConstantVolumeTriangleForce_fn(_system.Numfaces,
                                              &_system.vertices[0],
                                              &_system.faces[0],
                                              kappa_v,
                                              target_volume,
                                              current_volume,
                                              _system.get_box());
}