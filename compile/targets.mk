#
# License-Identifier: GPL
#
# Copyright (C) 2025 The Yambo Team
#
# Authors (see AUTHORS file for details): AM MN RR
#
# Build the standalone lelphc binary from  LetzElPhC C source.
# Requires yambo to have been configured first (config/setup must exist).
# Usage: make lelphc          -> builds plugins/ep/services/lelphc
#        make clean_lelphc    -> cleans the standalone build

lelphc:
	$(MAKE) -C $(srcdir)/plugins/ep/services

clean_lelphc:
	-$(MAKE) -C $(srcdir)/plugins/ep/services clean
