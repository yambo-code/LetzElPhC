// Funtions related to timings in the code
//
//

#pragma once

void init_ELPH_clocks(void);

void ELPH_start_clock(const char* str);

void ELPH_stop_clock(const char* str);

void cleanup_ELPH_clocks(void);

void print_ELPH_clock_summary(void);
