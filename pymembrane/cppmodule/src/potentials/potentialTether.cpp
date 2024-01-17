#include "potentialTether.hpp"

void ComputeVertexLimitForce_fn(int Numvertices,
                                HE_VertexProp *vertices,
                                real *_Fx,
                                real *_Fy,
                                real *_Fz)
{
    for (int vertex_index = 0; vertex_index < Numvertices; vertex_index++)
    {
        int type = vertices[vertex_index].type;

        vertices[vertex_index].forceC.x += _Fx[type];
        vertices[vertex_index].forceC.y += _Fy[type];
        vertices[vertex_index].forceC.z += _Fz[type];
    }
}
void ComputeVertexTetherEnergy::compute(void)
{

    ComputeVertexLimitForce_fn(_system.Numvertices,
                               &_system.vertices[0],
                               &m_Fx[0],
                               &m_Fy[0],
                               &m_Fz[0]);
}
