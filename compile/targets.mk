#
# License-Identifier: GPL
#
# Copyright (C) 2025 The Yambo Team
#
# Authors (see AUTHORS file for details): AM MN RR
#
lelphc:
	$(MAKE) -C $(srcdir)/plugins/letz/services

clean_lelphc:
	-$(MAKE) -C $(srcdir)/plugins/letz/services clean
