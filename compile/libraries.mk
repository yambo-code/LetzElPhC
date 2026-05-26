#
# License-Identifier: GPL
#
# Copyright (C) 2025 The Yambo Team
#
# Authors (see AUTHORS file for details): AM MN
#
ifeq ($(wildcard compile/global/defs.mk),compile/global/defs.mk)
  include compile/defs.mk
endif
#
# Library link order is INTRINSIC to LetzElPhC's dependency.
# Static linker is left-to-right: symbol-users before symbol-providers.
#
# Dependency tiers (each tier depends only on equal/lower tiers):
#   T5 drivers      : interpolation, elph
#   T4 QE I/O       : io/qe
#   T3 main I/O     : io
#   T2 kernels+pp   : dvloc, nonloc, symmetries, preprocessor
#   T1 primitives   : wfc, fft, phonon, parser
#   T0 base+vendors : common, common/cwalk, common/ELPH_hash_map, common/kdtree,
#                     io/ezxml, parser/inih
#
# elph (T5) hosts the Fortran-to-C ABI boundary (ep_f2c_bridge.c):
#   elph_driver_f2c / elph_driver_cb_f2c / elph_driver_cb2_f2c
# These are the ONLY entry points yambo_ep calls into LetzElPhC.
#
EP_S1= [Sletz,interpolation] [Sletz,elph]
EP_S2= [Sletz,io/qe]
EP_S3= [Sletz,io]
EP_S4= [Sletz,dvloc] [Sletz,nonloc] [Sletz,symmetries] [Sletz,preprocessor]
EP_S5= [Sletz,wfc] [Sletz,fft] [Sletz,phonon] [Sletz,parser]
EP_S6= [Sletz,common] [Sletz,common/cwalk] [Sletz,common/ELPH_hash_map] \
        [Sletz,common/kdtree] [Sletz,io/ezxml] [Sletz,parser/inih]
#
SERVICES_LIBS += ${EP_S1} ${EP_S2} ${EP_S3} ${EP_S4} ${EP_S5} ${EP_S6}
