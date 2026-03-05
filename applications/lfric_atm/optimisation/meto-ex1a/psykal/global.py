##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################


'''
PSyclone transformation script for the LFRic (Dynamo0p3) API to apply
colouring, OpenMP and redundant computation to the level-1 halo for
the initialisation built-ins generically.

'''

from psyclone_tools import (redundant_computation_setval, colour_loops,
                            openmp_parallelise_loops,
                            view_transformed_schedule)


def trans(psyir):
    '''
    Applies PSyclone colouring, OpenMP and redundant computation
    transformations.

    :param psyir: the PSyIR of the PSy-layer.
    :type psyir: :py:class:`psyclone.psyir.nodes.FileContainer`

    '''


    # Restrict tiling to kernels that can benefit
    tiling_kernels = [
        'matrix_vector_code',
        'pressure_gradient_p0_code',
        'nodal_xyz_coordinates_code',
        'schur_backsub_code',
        'apply_mixed_operator_code'#,
        #'apply_mixed_lu_operator_code',
        #'apply_elim_mixed_lp_operator_code',
        #'scaled_matrix_vector_code'
        #'opt_apply_variable_hx_code',
        #'dg_matrix_vector_code',
        #'dg_inc_matrix_vector_code',
        #'helmholtz_operator_code',
        #'apply_helmholtz_operator_code'#,
    #    'held_suarez_fv_wind_code'
    ]

    redundant_computation_setval(psyir)
    colour_loops(psyir,enable_tiling=True,tiling_kernel_list=None)
    openmp_parallelise_loops(psyir,enable_profiler=True)
    view_transformed_schedule(psyir)
