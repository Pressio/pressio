
#------------------------------------------------------------------------
# to include the lumped headers to have access to package classes
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../..)

# include the utilities in utils (adding this here means
# all children cmakelists will inherit this include )
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/utils)

ADD_SUBDIRECTORY(burgers1d)
ADD_SUBDIRECTORY(steady_linear_adv_diff_1d)
ADD_SUBDIRECTORY(steady_linear_adv_diff_2d)
ADD_SUBDIRECTORY(unsteady_linear_adv_diff_1d)
ADD_SUBDIRECTORY(unsteady_nonlinear_adv_diff_reaction_2d)
ADD_SUBDIRECTORY(unsteady_nonlinear_adv_diff_reaction_flame_2d)
