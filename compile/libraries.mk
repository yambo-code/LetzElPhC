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
SERVICES_LIBS += [Sep,common/kdtree] [Sep,phonon] [Sep,common/ELPH_hash_map] [Sep,interpolation] 
SERVICES_LIBS +=  [Sep,common/cwalk] [Sep,fft] [Sep,symmetries] [Sep,wfc] 
SERVICES_LIBS += [Sep,common] [Sep,preprocessor] [Sep,dvloc] [Sep,parser/inih] [Sep,parser] [Sep,io/ezxml] [Sep,io/qe] [Sep,nonloc] [Sep,io] [Sep,elph]
#
