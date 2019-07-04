dnl -*- autoconf -*-
# Macros needed to find transport and dependent libraries.  They are called by
# the macros in ${top_src_dir}/dependencies.m4, which is generated by
# "dunecontrol autogen"

# Additional checks needed to build transport
# This macro should be invoked by every module which depends on transport, as
# well as by transport itself
AC_DEFUN([TRANSPORT_CHECKS])

# Additional checks needed to find transport
# This macro should be invoked by every module which depends on transport, but
# not by transport itself
AC_DEFUN([TRANSPORT_CHECK_MODULE],
[
  DUNE_CHECK_MODULES([transport],[transport/transport.hh])
])
