#
# License-Identifier: GPL
#
# Copyright (C) 2016 The Yambo Team
#
# Authors (see AUTHORS file for details): HM
#
# VARIABLES (static and dynamical)
#
include dyn_variables.mk
include $(compdir)/compile/local/static_variables.mk
include $(compdir)/compile/local/defs.mk
include $(compdir)/compile/defs.mk
#
# SETUP
#
include $(compdir)/config/setup
#
# Module dependencies
#
-include local_modules.dep
#
# Include directories (headers)
#
include $(compdir)/compile/headers.mk
#
# Libraries 
#
include $(compdir)/compile/libraries.mk
ifneq "$(wildcard $(compdir)/plugins/plugins.pulled)" ""
 include $(compdir)/plugins/*/compile/headers_and_libs.mk
endif
#
# OBJECTS
#
include objects.mk
#
# OPERATIONS
#
include $(compdir)/compile/local/operations.mk
#
# RULES
#
include $(compdir)/compile/local/rules.mk
#
# FUNCTIONS
#
include $(compdir)/compile/local/functions.mk
#
