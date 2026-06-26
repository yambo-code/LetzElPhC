#
# License-Identifier: GPL
#
# Copyright (C) 2025 The Yambo Team
#
# Authors (see AUTHORS file for details): RR AM
#
# Help blurb for the ep (LetzElPhC) plugin, called from
# compile/global/functions/help.mk when plugins/ep.pulled is present.
#
define ep_help
 $(ECHO)  " ep   plugin =  LetzElPhC electron-phonon coupling plugin";\
 $(ECHO)  "                Trigger: yambo -collisions ep";\
 $(ECHO)  "                Source: plugins/ep/services/, integration: packages/el-ph/";\
 $(ECHO)
endef
