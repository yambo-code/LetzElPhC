#
# License-Identifier: GPL
#
# Copyright (C) 2025 The Yambo Team
#
# Authors (see AUTHORS file for details): RR AM
#
# LetzElPhC ep plugin - extra include paths for the vendored C source.
# C files under plugins/ep/services/<sub>/*.c use relative includes
# like "common/dtypes.h", so the plugin services root must be on the
# include path.
#
iheaders        += $(IFLAG)$(srcdir)/plugins/letz/services
include_headers += $(IFLAG)$(srcdir)/plugins/letz/services
llibs           += -lcblas
