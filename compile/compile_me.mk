#
# License-Identifier: GPL
#
# Copyright (C) 2025 The Yambo Team
#
# Authors (see AUTHORS file for details): MN
#
# Co-Authors (see AUTHORS file for details): RR AM
# LetzElPhC ep plugin - yambo build system stub
# Integration handled via SERVICES_LIBS in libraries.mk
#
ifneq (,$(findstring yambo_ep,$(MAKECMDGOALS)))
 Y_PRECMP= -D_Y6_LETZ
endif
