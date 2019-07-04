add_executable("pdelab_ss" simple.cc driver_simple.hh operator_simple.hh bctype_simple.hh)
target_link_dune_default_libraries("pdelab_ss")
