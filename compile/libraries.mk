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
EP_S1= [Sletz,nonloc] [Sletz,io/ezxml] [Sletz,parser/inih] [Sletz,io] 
EP_S2= [Sletz,io/qe] [Sletz,common] [Sletz,dvloc] [Sletz,elph] [Sletz,wfc] [Sletz,symmetries] [Sletz,fft]
EP_S3= [Sletz,common/cwalk] [Sletz,preprocessor] [Sletz,interpolation] [Sletz,common/ELPH_hash_map] [Sletz,phonon] [Sletz,common/kdtree] [Sletz,parser]
#
ifneq (,$(filter yambo_ep ypp_ep,$(MAKECMDGOALS)))
  SERVICES_LIBS += ${EP_S1} ${EP_S2} ${EP_S3}
endif
