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
EP_S1= [Sep,nonloc] [Sep,io/ezxml] [Sep,parser/inih] [Sep,io] 
EP_S2= [Sep,io/qe] [Sep,common] [Sep,dvloc] [Sep,elph] [Sep,wfc] [Sep,symmetries] [Sep,fft]
EP_S3= [Sep,common/cwalk] [Sep,preprocessor] [Sep,interpolation] [Sep,common/ELPH_hash_map] [Sep,phonon] [Sep,common/kdtree] [Sep,parser]
#
SERVICES_LIBS += ${EP_S1} ${EP_S2} ${EP_S3}

PLUGINS_include = test
