#include "potentialConstantAreaGlobal.hpp"

void ComputeVertexConstantGlobalAreaEnergy::recalculate_mesh_area(void)
{
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

real ComputeVertexConstantGlobalAreaEnergy_lambda(const real _current_area,
                                                    const real _kappa_ag,
                                                    const real _target_area)
{
    real energy = _kappa_ag/_target_area * (_current_area - _target_area) * (_current_area - _target_area);
    return energy;
}

void ComputeVertexConstantGlobalAreaEnergy_fn(const int Numfaces,
                                                HE_FaceProp *faces,
                                                const real _current_area,
                                                const real _kappa_ag,
                                                const real _target_area)
{

    real energy = ComputeVertexConstantGlobalAreaEnergy_lambda(_current_area, _kappa_ag, _target_area);
    for (int face_index = 0; face_index < Numfaces; face_index++)
    {
        /// Add energy to that face
        faces[face_index].energy += energy/Numfaces;
    }
}

void ComputeVertexConstantGlobalAreaEnergy::compute_energy(void)
{

    ComputeVertexConstantGlobalAreaEnergy_fn(_system.Numfaces,
                                               &_system.faces[0],
                                               current_area,
                                               kappa_ag,
                                               target_area);
}

real ComputeVertexConstantGlobalAreaEnergy::compute_edge_energy(int query_edge_index)
{
    // reset area to recalculate
    this->recalculate_mesh_area();
    // we need to loop the 2 faces that are connected to the edge_index
    pymemb::vector<int> face_vec{_system.edges[query_edge_index].face_k, _system.edges[query_edge_index].face_l};
    // reset energy
    real edge_energy = 0.0;
    for (auto face_index : face_vec)
    {
        edge_energy += ComputeVertexConstantGlobalAreaEnergy_lambda(current_area, kappa_ag, target_area);
    }
    return edge_energy;
}

real ComputeVertexConstantGlobalAreaEnergy::compute_vertex_energy(int query_vertex_index)
{
    // reset area to recalculate
    this->recalculate_mesh_area();
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
            energy += kappa_ag*(current_area-target_area)*(current_area-target_area)/(target_area);
        }

        int he_prev = _system.halfedges[he].prev;
        he = _system.halfedges[he_prev].pair;
    } while ((he != first));
    return energy;
}

void ComputeVertexConstantAreaTriangleForce_fn(const int Numfaces,
                                               HE_VertexProp *vertices,
                                               const HE_FaceProp *faces,
                                               const real _kappa_ag,
                                               const real _target_area,
                                               const real _current_area,
                                               const BoxType _box)
{

    real3 r1;
    real3 r2;
    real3 r3;
    real3 cross;

    int v1, v2, v3;

    real force_factor = -_kappa_ag*2*(_current_area-_target_area)/(_target_area);
    
    for (int face_index = 0; face_index < Numfaces; face_index++)
    {
        int type = faces[face_index].type;
        v1 = faces[face_index].v1;
        v2 = faces[face_index].v2;
        v3 = faces[face_index].v3;

        r1 = vertices[v1].r;
        r2 = vertices[v2].r;
        r3 = vertices[v3].r;

        auto r12 = pymemb::vector_subtract(r2, r1, _box);
        auto r13 = pymemb::vector_subtract(r3, r1, _box);
        auto r32 = pymemb::vector_subtract(r2, r3, _box);

        // we use dA/dvertex = 1/2 normal x u, see www.cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf Appendix A
        // The surface should have been set up so that this normal is oriented away from the volume.
        real3 normal = pymemb::compute_normal_triangle_unit(r1, r2, r3, _box);

        /*----------------------------------------------------------------------------------------------------------------*/
        /*-----------------------------------           ACTUAL CALCULATION        ----------------------------------------*/
        /*----------------------------------------------------------------------------------------------------------------*/
        // v1
        cross = pymemb::vector_cross(r32, normal);
        vertices[v1].forceC.x += 0.5*force_factor * cross.x;
        vertices[v1].forceC.y += 0.5*force_factor * cross.y;
        vertices[v1].forceC.z += 0.5*force_factor * cross.z;

        // v2
        cross = pymemb::vector_cross(r13, normal);
        vertices[v2].forceC.x += 0.5*force_factor * cross.x;
        vertices[v2].forceC.y += 0.5*force_factor * cross.y;
        vertices[v2].forceC.z += 0.5*force_factor * cross.z;

        // v3, be careful of direction!
        cross = pymemb::vector_cross(normal, r12);
        vertices[v3].forceC.x += 0.5*force_factor * cross.x;
        vertices[v3].forceC.y += 0.5*force_factor * cross.y;
        vertices[v3].forceC.z += 0.5*force_factor * cross.z;
    }
}

void ComputeVertexConstantGlobalAreaEnergy::compute(void)
{

    this->recalculate_mesh_area();

    ComputeVertexConstantAreaTriangleForce_fn(_system.Numfaces,
                                              &_system.vertices[0],
                                              &_system.faces[0],
                                              kappa_ag,
                                              target_area,
                                              current_area,
                                              _system.get_box());
}